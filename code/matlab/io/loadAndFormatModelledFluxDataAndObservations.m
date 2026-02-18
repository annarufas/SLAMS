function [modFluxMonthlyProfile,modFluxMonthlyDh,modFluxMonthlyDhDepths,...
          modFluxSeasonalProfile,modFluxSeasonalDh,modFluxSeasonalDhDepths,... 
          obsDataMonthlyNoOutliers,obsErrorMonthlyNoOutliers,obsDataSeasonal,...
          obsErrorSeasonal,localBicepExFluxMonthly,bicepExDepthMonthly,...
          localBicepExFluxSeasonal,bicepExDepthSeasonal] =...
    loadAndFormatModelledFluxDataAndObservations(choiceTypeGridDomain,config,...
        fullpathModelRunsDir,filenameSlamsOutput,filenameObsPocFlux,...
        filenameObsPicAndBsiFlux,filenameTimeseriesInformation,...
        filenameBicepExportFluxDunne,filenameBicepExportFluxHenson,...
        filenameBicepExportFluxLi,filenameZeu)

% ======================================================================= %
%                                                                         %
% Loads and formats modelled and observational local monthly flux data in %
% preparation for plotting or model optimisation.                         %
%                                                                         %
% ======================================================================= %

idxSeason = config.seasonIndices;

%% SECTION 1: Load time-series metadata from observational compilation

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
nLocs = length(STATION_NAMES);

% Rearrange model locations to match observational order
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'}; % as in config.gridLats, config.gridLons
[~,reorderModLocIdx] = ismember(STATION_NAMES,currentModLocationOrder); % get reordering indices

% Log progress
logID = fopen(fullfile('.','data','interim','log_zeu_calc.log'),'w'); 

%% SECTION 2: Load observational data from compilation

% POC flux compilation 
load(fullfile('.','data','processed',filenameObsPocFlux),...
    'pocFluxMonthlyDhAvg','pocFluxMonthlyDhErrTot','LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED')

% PIC and bSi flux compilation
load(fullfile('.','data','processed',filenameObsPicAndBsiFlux),...
    'picbsiFluxMonthlyDhAvg','picbsiFluxMonthlyDhErrTot','LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED')

nTargetDepths = size(LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED,4);

% Transform units of observed data to match model units:
% mmol POC m-2 d-1 --> mg C m-2 d-1 
% mmol PIC m-2 d-1 --> mg CaCO3 m-2d-1
% mmol Si m-2 d-1 --> mg opal m-2 d-1 
obsDataMonthly = {...
    config.molarMassCarbon.*pocFluxMonthlyDhAvg,...
    config.molarMassCaCO3.*picbsiFluxMonthlyDhAvg(:,:,:,1),...
    config.molarMassBiogenicSilica.*picbsiFluxMonthlyDhAvg(:,:,:,2)
};
obsErrorMonthly = {...
    config.molarMassCarbon.*pocFluxMonthlyDhErrTot,...
    config.molarMassCaCO3.*picbsiFluxMonthlyDhErrTot(:,:,:,1),...
    config.molarMassBiogenicSilica.*picbsiFluxMonthlyDhErrTot(:,:,:,2)
};

% Outlier removal in observations
[obsDataMonthlyNoOutliers,obsErrorMonthlyNoOutliers] = removeOutliersInFluxObservations(...
    obsDataMonthly,obsErrorMonthly,nTargetDepths,nLocs);

% Calculate seasonal average
obsDataSeasonal = cell(1,3);
obsErrorSeasonal = cell(1,3);

for iTracer = 1:3
    obsDataTracer = obsDataMonthlyNoOutliers{iTracer}; % nTargetDepths x 12 months x nLocs
    obsErrorTracer = obsErrorMonthlyNoOutliers{iTracer}; 

    % Preallocate seasonal arrays
    dataSeasonal  = NaN(nTargetDepths,4,nLocs);
    errorSeasonal = NaN(nTargetDepths,4,nLocs);

    % Loop through seasons
    for iSeason = 1:4
        dataSeasonal(:,iSeason,:)  = mean(obsDataTracer(:,idxSeason{iSeason},:),2,'omitnan');
        errorSeasonal(:,iSeason,:) = mean(obsErrorTracer(:,idxSeason{iSeason},:),2,'omitnan');
    end

    % Store seasonal results
    obsDataSeasonal{iTracer}  = dataSeasonal;
    obsErrorSeasonal{iTracer} = errorSeasonal;
end 

%% SECTION 3: Load observational data from BICEP export flux climatology

load(fullfile('.','data','interim',filenameBicepExportFluxDunne),'ef_dunne_avg','ef_lon','ef_lat')
load(fullfile('.','data','interim',filenameBicepExportFluxHenson),'ef_henson_avg')
load(fullfile('.','data','interim',filenameBicepExportFluxLi),'ef_li_avg')
globalBicepExFluxMonthly = NaN([size(ef_dunne_avg),3]);
globalBicepExFluxMonthly(:,:,:,1) = ef_dunne_avg;
globalBicepExFluxMonthly(:,:,:,2) = ef_henson_avg;
globalBicepExFluxMonthly(:,:,:,3) = ef_li_avg;

% Extract data at our locations of interest
localBicepExFluxMonthly = NaN(config.nLocs,12,3); % mg C m-2 d-1
for iLoc = 1:config.nLocs
    [~,iLon] = min(abs(ef_lon(:) - config.lons(iLoc)));
    [~,iLat] = min(abs(ef_lat(:) - config.lats(iLoc)));
    localBicepExFluxMonthly(iLoc,:,:) = globalBicepExFluxMonthly(iLat,iLon,:,:);
end  

% Rearrange model locations to match compilation observational order
localBicepExFluxMonthly = localBicepExFluxMonthly(reorderModLocIdx,:,:);

% Calculate seasonal average
localBicepExFluxSeasonal = NaN(nLocs,4,3);
for iSeason = 1:4
    localBicepExFluxSeasonal(:,iSeason,:) = mean(...
        localBicepExFluxMonthly(:,idxSeason{iSeason},:),2,'omitnan');
end

%% SECTION 4: Get euphotic layer depth (i.e., export depth) for the BICEP product

% Load the global-ocean euphotic layer depth product
load(fullfile('.','data','interim',filenameZeu),'zeu','zeu_lat','zeu_lon')
 
% Query points for interpolation
qLats = config.lats;
qLons = config.lons;

% Original data grid
[X,Y,T] = ndgrid(zeu_lat,zeu_lon,(1:12)');

% Interpolant 
F = griddedInterpolant(X, Y, T, zeu, 'linear'); 

% Extract data for the study locations (defined by qLats and qLons)
qZeuMonthly = NaN(12,config.nLocs);
for iLoc = 1:config.nLocs
    [qX,qY,qT] = ndgrid(qLats(iLoc),qLons(iLoc),(1:12)');

    % Perform interpolation of the dataset to the current location's coordinates
    localZeu = F(qX,qY,qT);

    % Replace NaNs using an interpolation method
    if any(isnan(localZeu)) && sum(~isnan(localZeu)) > 2 
        localZeuFilled = fillNansStrategicallyInSurfaceOceanData(localZeu,(1:12),logID);
        localZeuFilled(isnan(localZeuFilled)) = 0;
        qZeuMonthly(:,iLoc) = localZeuFilled;
    else
        qZeuMonthly(:,iLoc) = localZeu;
    end
end

% Rearrange model locations to match compilation observational order
bicepExDepthMonthly = qZeuMonthly(:,reorderModLocIdx);

% Calculate seasonal values
bicepExDepthSeasonal = NaN(4,nLocs);
for iSeason = 1:4
    bicepExDepthSeasonal(iSeason,:) = mean(...
        bicepExDepthMonthly(idxSeason{iSeason},:),1,'omitnan');
end

%% SECTION 5: Load modelled data

load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')

% Extract flux components
tracers = struct(...
    'POC', squeeze(output.flux.monthly(1,:,:,:)) + squeeze(output.flux.monthly(2,:,:,:)),...
    'PIC', squeeze(output.flux.monthly(3,:,:,:)),...
    'BSi', squeeze(output.flux.monthly(4,:,:,:))...
);
tracerNames = fieldnames(tracers);

% Assemble modelled flux array (nDepths x 12 months x nLocs x 3 tracers)
modFluxMonthlyProfile = NaN([size(tracers.POC),numel(tracerNames)]); 
for iTracer = 1:numel(tracerNames)
    modFluxMonthlyProfile(:,:,:,iTracer) = tracers.(tracerNames{iTracer}); 
end 
modFluxMonthlyProfile = modFluxMonthlyProfile(:,:,reorderModLocIdx,:);

% Calculate seasonal average
modFluxSeasonalProfile = NaN(size(modFluxMonthlyProfile,1),4,nLocs,numel(tracerNames));
for iSeason = 1:4
    modFluxSeasonalProfile(:,iSeason,:,:) = mean(modFluxMonthlyProfile(:,idxSeason{iSeason},:,:),2,'omitnan');
end

%% SECTION 6: Extract modelled data at observation depth horizons

obsFluxDepthsList = cat(5, ...
    LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED, ...
    LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(:,:,:,:,1), ...
    LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(:,:,:,:,2));

% Reorder dimensions before passing in to the next functions 
% (12 x 6 x 2 x 4 x 3 --> 4 x 12 x 6 x 2 x 3)
obsFluxDepthsList = permute(obsFluxDepthsList, [4, 1, 2, 3, 5]);

[modFluxMonthlyDh,modFluxMonthlyDhDepths,~,~] =...
    extractKeyDepthHorizonDataInModelledFluxes(config,modFluxMonthlyProfile,obsFluxDepthsList);

% Calculate seasonal average
modFluxSeasonalDh = NaN(nTargetDepths,4,nLocs,2,numel(tracerNames));
modFluxSeasonalDhDepths = NaN(nTargetDepths,4,nLocs,numel(tracerNames));
for iSeason = 1:4
    modFluxSeasonalDh(:,iSeason,:,:,:) = mean(modFluxMonthlyDh(:,idxSeason{iSeason},:,:,:),2,'omitnan');
    modFluxSeasonalDhDepths(:,iSeason,:,:) = mean(modFluxMonthlyDhDepths(:,idxSeason{iSeason},:,:),2,'omitnan');
end

end % loadAndFormatModelledFluxDataAndObservations