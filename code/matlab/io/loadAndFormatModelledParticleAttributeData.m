function [particleAttributesByScAnnualDhMetadata,...
          particleAttributesByScSeasonalDhMetadata,...
          particleAttributesByTypeAnnualDhMetadata,...
          fluxMonthlyDh,fluxSeasonalDh] =...
    loadAndFormatModelledParticleAttributeData(config,filenameTimeseriesInformation,...
        filenameObsPocFlux,fullpathModelRunsDir,filenameSlamsOutput)

idxSeason = config.seasonIndices;

%% SECTION 1: Load model output data

load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')

%% SECTION 2: Load metadata

% Load station-related information from POC flux observation compilation
% (this way, the order of locations will match that of plotted flux data 
% and particle number data)
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};

% Rearrange model locations to match observational order
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'}; % as in config.gridLats, config.gridLons
[~,reorderModLocIdx] = ismember(desiredLocationOrder,currentModLocationOrder); % get reordering indices

%% SECTION 3: Depth selection for plotting

% POC flux compilation depths
load(fullfile('.','data','processed',filenameObsPocFlux),'LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED')

% Rearrange observed locations depths
[~,reorderObsLocIdx] = ismember(desiredLocationOrder,STATION_NAMES); % get reordering indices
LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(:,reorderObsLocIdx,:,:);

% Calculate annual average
obsDepthsAnnualAvg = squeeze(mean(LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED,1,'omitnan'));
obsDepthsMonthlyAvg = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED;

% Calculate seasonal average
obsDepthsSeasonalAvg = NaN(4,size(obsDepthsAnnualAvg,1),size(obsDepthsAnnualAvg,2),size(obsDepthsAnnualAvg,3));
for iSeason = 1:4
    obsDepthsSeasonalAvg(iSeason,:,:,:) = mean(obsDepthsMonthlyAvg(idxSeason{iSeason},:,:,:),1,'omitnan');
end

% Replace the last observed depth in sediment traps with a depth more
% relevant for particle number observations (last non-NaN)
localImagedDepths = output.deployment.imagSys;
for iLoc = 1:config.nLocs
    idxLastObsDepth = find(~isnan(localImagedDepths(:,iLoc)), 1, 'last'); % last non-NaN
    obsDepthsAnnualAvg(iLoc,:,4) = localImagedDepths(idxLastObsDepth,iLoc);
    obsDepthsMonthlyAvg(:,iLoc,:,4) = localImagedDepths(idxLastObsDepth,iLoc);
    obsDepthsSeasonalAvg(:,iLoc,:,4) = localImagedDepths(idxLastObsDepth,iLoc);
end

%% SECTION 4: Extract modelled particle attribute data by main size class at key depth horizons

% Extract ALL attribute data and reorder
particleAvgAttsAnnualProfile = output.particle.avgAttScAtDepthAnnual; % 3 sizes x nDepths x nAtts x nLocs
particleAvgAttsAnnualProfile = particleAvgAttsAnnualProfile(:,:,:,reorderModLocIdx);

particleAvgAttsSeasonalProfile = output.seasonal.particle.avgAttScAtDepth; % 3 sizes x nDepths x nAtts x 4 seasons x nLocs
particleAvgAttsSeasonalProfile = particleAvgAttsSeasonalProfile(:,:,:,:,reorderModLocIdx);

% Extract particle attribute data at depths of interest
nTargetDepths = size(obsDepthsAnnualAvg,3);
nImagedDepths = size(particleAvgAttsAnnualProfile,2);
nPartAtts = size(particleAvgAttsAnnualProfile,3);
nPartSizes = size(particleAvgAttsAnnualProfile,1);
particleAvgAttsAnnualDh = NaN(nPartSizes,nTargetDepths,nPartAtts,config.nLocs);
particleAvgAttsSeasonalDh = NaN(nPartSizes,nTargetDepths,nPartAtts,4,config.nLocs);

validDepthMaskAnnual = false(nPartSizes,nImagedDepths,nPartAtts,config.nLocs,nTargetDepths);
validDepthMaskSeasonal = false(nPartSizes,nImagedDepths,nPartAtts,4,config.nLocs,nTargetDepths);
for iDh = 1:nTargetDepths
    for iLoc = 1:config.nLocs
        depthRangeTop = obsDepthsAnnualAvg(iLoc,1,iDh);    % upper boundary depth
        depthRangeBottom = obsDepthsAnnualAvg(iLoc,2,iDh); % lower boundary depth

        % Only proceed if both are numbers rather than NaN
        if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
            % Find the closest depth greater than or equal to depthRangeTop
            [~,idxTop] = min(abs(localImagedDepths(:,iLoc) - depthRangeTop));

            % Find the closest depth less than or equal to depthRangeBottom
            [~,idxBottom] = min(abs(localImagedDepths(:,iLoc) - depthRangeBottom));
            if idxBottom > nImagedDepths
                idxBottom = nImagedDepths;
            end

            % Extract valid modelled depth range
            validDepthMaskAnnual(:,(idxTop:idxBottom),:,iLoc,iDh) = true;
            validDepthMaskSeasonal(:,(idxTop:idxBottom),:,:,iLoc,iDh) = true;
        end
    end % iLoc
end % iDh

% Extract particle attribute data
for iDh = 1:nTargetDepths 
    % Annual
    currentMaskAnnual = squeeze(validDepthMaskAnnual(:,:,:,:,iDh));
    partAttValuesAnnual = particleAvgAttsAnnualProfile; 
    partAttValuesAnnual(~currentMaskAnnual) = NaN;
    particleAvgAttsAnnualDh(:,iDh,:,:) = squeeze(mean(partAttValuesAnnual,2,'omitnan'));
    % Seasonal
    currentMaskSeasonal = squeeze(validDepthMaskSeasonal(:,:,:,:,:,iDh));
    partAttValuesSeasonal = particleAvgAttsSeasonalProfile; 
    partAttValuesSeasonal(~currentMaskSeasonal) = NaN;
    particleAvgAttsSeasonalDh(:,iDh,:,:,:) = squeeze(mean(partAttValuesSeasonal,2,'omitnan'));
end

%% SECTION 5: Extract modelled particle attribute data by particle type at key depth horizons

% Extract ALL attribute data and reorder
particleAvgAttInTyProfile = output.particle.avgAttInTyAnnual; % 8 types x nDepths x nAtts x nLocs
particleAvgAttInTyProfile = particleAvgAttInTyProfile(:,:,:,reorderModLocIdx);

% Extract particle attribute data at depths of interest
nPartTypes = size(particleAvgAttInTyProfile,1);
nPartAtts = size(particleAvgAttInTyProfile,3);
particleAvgAttInTyAnnualDh = NaN(nPartTypes,nTargetDepths,nPartAtts,config.nLocs);
particleAvgAttInTyAnnualDhDepths = NaN(size(particleAvgAttInTyAnnualDh));

validDepthMaskAnnual = false(nPartTypes,nImagedDepths,nPartAtts,config.nLocs,nTargetDepths);
for iDh = 1:nTargetDepths
    for iLoc = 1:config.nLocs
        depthRangeTop = obsDepthsAnnualAvg(iLoc,1,iDh);    % upper boundary depth
        depthRangeBottom = obsDepthsAnnualAvg(iLoc,2,iDh); % lower boundary depth

        % Only proceed if both are numbers rather than NaN
        if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
            % Find the closest depth greater than or equal to depthRangeTop
            [~,idxTop] = min(abs(localImagedDepths(:,iLoc) - depthRangeTop));

            % Find the closest depth less than or equal to depthRangeBottom
            [~,idxBottom] = min(abs(localImagedDepths(:,iLoc) - depthRangeBottom));
            if idxBottom > nImagedDepths
                idxBottom = nImagedDepths;
            end

            % Extract valid modelled depth range
            validDepthMaskAnnual(:,(idxTop:idxBottom),:,iLoc,iDh) = true;

            % Calculate the average monthly unique modelled depth
            particleAvgAttInTyAnnualDhDepths(:,iDh,:,iLoc) =...
                (localImagedDepths(idxTop,iLoc)+localImagedDepths(idxBottom,iLoc))/2;
        end
    end % iLoc
end % iDh

% Extract particle attribute data
for iDh = 1:nTargetDepths
    currentMaskAnnual = squeeze(validDepthMaskAnnual(:,:,:,:,iDh));
    partAttValuesAnnual = particleAvgAttInTyProfile(:,:,:,:); 
    partAttValuesAnnual(~currentMaskAnnual) = NaN;
    particleAvgAttInTyAnnualDh(:,iDh,:,:) = squeeze(mean(partAttValuesAnnual,2,'omitnan'));
end

%% SECTION 6: Extract modelled flux components at key depth horizons

% Extract flux components
tracers = struct(...
    'POC', squeeze(output.flux.monthly(1,:,:,:)) + squeeze(output.flux.monthly(2,:,:,:)),...
    'PIC', squeeze(output.flux.monthly(3,:,:,:)),...
    'BSi', squeeze(output.flux.monthly(4,:,:,:))...
);
tracerNames = fieldnames(tracers);

% Assemble modelled flux array (nDepths x 12 x nLocs x 3 tracers)
fluxMonthlyProfile = NaN([size(tracers.POC),numel(tracerNames)]); 
for iTracer = 1:numel(tracerNames)
    fluxMonthlyProfile(:,:,:,iTracer) = tracers.(tracerNames{iTracer}); 
end 
fluxMonthlyProfile = fluxMonthlyProfile(:,:,reorderModLocIdx,:);

% Extract flux data at depths of interest
localSedTrapDepths = output.deployment.sedTrap;
nTargetDepths = size(obsDepthsMonthlyAvg,4);
nSedTrapDepths = size(fluxMonthlyProfile,1);
fluxMonthlyDh = NaN(nTargetDepths,12,config.nLocs,numel(tracerNames));
fluxMonthlyDhDepths = NaN(size(fluxMonthlyDh));

validDepthMaskMonthly = false(nSedTrapDepths,12,config.nLocs,numel(tracerNames),nTargetDepths);
for iDh = 1:nTargetDepths
    for iMonth = 1:12
        for iLoc = 1:config.nLocs
            depthRangeTop = obsDepthsMonthlyAvg(iMonth,iLoc,1,iDh);    % upper boundary depth
            depthRangeBottom = obsDepthsMonthlyAvg(iMonth,iLoc,2,iDh); % lower boundary depth
    
            % Only proceed if both are numbers rather than NaN
            if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
                % Find the closest depth greater than or equal to depthRangeTop
                [~,idxTop] = min(abs(localSedTrapDepths(:,iLoc) - depthRangeTop));
    
                % Find the closest depth less than or equal to depthRangeBottom
                [~,idxBottom] = min(abs(localSedTrapDepths(:,iLoc) - depthRangeBottom));
                if idxBottom > nSedTrapDepths
                    idxBottom = nSedTrapDepths;
                end
    
                % Extract valid modelled depth range
                validDepthMaskMonthly((idxTop:idxBottom),iMonth,iLoc,:,iDh) = true;
    
                % Calculate the average monthly unique modelled depth
                fluxMonthlyDhDepths(iDh,iMonth,iLoc,:) =...
                    (localSedTrapDepths(idxTop,iLoc)+localSedTrapDepths(idxBottom,iLoc))/2;
            end
        end % iLoc
    end % iMonth
end % iDh

% Extract flux data
for iDh = 1:nTargetDepths
    currentMask = squeeze(validDepthMaskMonthly(:,:,:,:,iDh));
    fluxValues = fluxMonthlyProfile; 
    fluxValues(~currentMask) = NaN;
    fluxMonthlyDh(iDh,:,:,:) = squeeze(mean(fluxValues,1,'omitnan'));
end

% Calculate seasonal average
fluxSeasonalDh = NaN(size(fluxMonthlyDh,1),4,config.nLocs,numel(tracerNames));
for iSeason = 1:4
    fluxSeasonalDh(:,iSeason,:,:) = mean(fluxMonthlyDh(:,idxSeason{iSeason},:,:),2,'omitnan');
end

%% SECTION 7: Assemble particle attribute into a structure for easier management of the array

% Got rid of attribute #9, depth
particleAttributesByScAnnualDhMetadata = {... % nLocs x 2 sizes x 4 depths
    'poc',   '$${\mathrm{POC}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,10,:)), [3,1,2]), 'POC content (pmol)';
    'tepc',  '$${\mathrm{TEP-C}}$$',            permute(squeeze(particleAvgAttsAnnualDh(2:3,:,11,:)), [3,1,2]), 'TEP-C content (pmol)';
    'alpha', '$$\alpha_{\mathrm{p}}$$',         permute(squeeze(particleAvgAttsAnnualDh(2:3,:,4,:)), [3,1,2]),  'Stickiness';
    'calc',  '$${\mathrm{CaCO_{3}}}$$',         permute(squeeze(particleAvgAttsAnnualDh(2:3,:,13,:)), [3,1,2]), 'CaCO3 content (pmol)';
    'bsi',   '$${\mathrm{bSi}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,12,:)), [3,1,2]), 'Opal content (pmol)';
    'clay',  '$${\mathrm{Clay}}$$',             permute(squeeze(particleAvgAttsAnnualDh(2:3,:,14,:)), [3,1,2]), 'Clay content (pmol)';    
    'vol',   '$$V_{\mathrm{p}}^{\mathrm{S}}$$', permute(squeeze(particleAvgAttsAnnualDh(2:3,:,6,:)), [3,1,2]),  'Solid volume (\mu^{3})';
    'd3',    '$$D_{\mathrm{3}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,7,:)), [3,1,2]),  'Fractal dimension';
    'esd',   '$${\mathrm{ESD}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,15,:)), [3,1,2]), 'ESD (\mum)';
    'rpp',   '$${r_{\mathrm{pp}}$$',            permute(squeeze(particleAvgAttsAnnualDh(2:3,:,8,:)), [3,1,2]),  'Radius primary particle (\mum)';
    'por',   '$$P_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,5,:)), [3,1,2]),  'Porosity';
    'rho',   '$$\rho_{\mathrm{p}}$$',           permute(squeeze(particleAvgAttsAnnualDh(2:3,:,2,:)), [3,1,2]),  'Density (g cm^{-3})';
    'vsink', '$$v_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,3,:)), [3,1,2]),  'Sinking velocity (m d^{-1})';
    'pnum',  '$$N$$',                           permute(squeeze(particleAvgAttsAnnualDh(2:3,:,1,:)), [3,1,2]),  'Particle number (# L^{-1})';
};

particleAttributesByScSeasonalDhMetadata = {... % nLocs x 2 sizes x 4 depths x 4 seasons
    'poc',   '$${\mathrm{POC}}$$',              permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,10,:,:)), [4,1,2,3]), 'POC content (pmol)';
    'tepc',  '$${\mathrm{TEP-C}}$$',            permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,11,:,:)), [4,1,2,3]), 'TEP-C content (pmol)';
    'alpha', '$$\alpha_{\mathrm{p}}$$',         permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,4,:,:)), [4,1,2,3]),  'Stickiness';
    'calc',  '$${\mathrm{CaCO_{3}}}$$',         permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,13,:,:)), [4,1,2,3]), 'CaCO3 content (pmol)';
    'bsi',   '$${\mathrm{bSi}}$$',              permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,12,:,:)), [4,1,2,3]), 'Opal content (pmol)';
    'clay',  '$${\mathrm{Clay}}$$',             permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,14,:,:)), [4,1,2,3]), 'Clay content (pmol)';    
    'vol',   '$$V_{\mathrm{p}}^{\mathrm{S}}$$', permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,6,:,:)), [4,1,2,3]),  'Solid volume (\mu^{3})';
    'd3',    '$$D_{\mathrm{3}}$$',              permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,7,:,:)), [4,1,2,3]),  'Fractal dimension';
    'esd',   '$${\mathrm{ESD}}$$',              permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,15,:,:)), [4,1,2,3]), 'ESD (\mum)';
    'rpp',   '$${r_{\mathrm{pp}}$$',            permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,8,:,:)), [4,1,2,3]),  'Radius primary particle (\mum)';
    'por',   '$$P_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,5,:,:)), [4,1,2,3]),  'Porosity';
    'rho',   '$$\rho_{\mathrm{p}}$$',           permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,2,:,:)), [4,1,2,3]),  'Density (g cm^{-3})';
    'vsink', '$$v_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,3,:,:)), [4,1,2,3]),  'Sinking velocity (m d^{-1})';
    'pnum',  '$$N$$',                           permute(squeeze(particleAvgAttsSeasonalDh(2:3,:,1,:,:)), [4,1,2,3]),  'Particle number (# L^{-1})';
};

particleAttributesByTypeAnnualDhMetadata = {... % nLocs x 8 types x 4 depths
    'esd',   '$${\mathrm{ESD}}$$',              permute(squeeze(particleAvgAttInTyAnnualDh(:,:,2,:)), [3,1,2]), 'ESD (\mum)';
    'por',   '$$P_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttInTyAnnualDh(:,:,5,:)), [3,1,2]), 'Porosity';
    'rho',   '$$\rho_{\mathrm{p}}$$',           permute(squeeze(particleAvgAttInTyAnnualDh(:,:,3,:)), [3,1,2]), 'Density (g cm^{-3})';
    'vsink', '$$v_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttInTyAnnualDh(:,:,4,:)), [3,1,2]), 'Sinking velocity (m d^{-1})';
};

end % loadAndFormatModelledParticleAttributeData