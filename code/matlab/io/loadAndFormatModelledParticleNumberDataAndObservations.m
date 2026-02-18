function [modPnumByScMonthlyProfile,modPnumByScSeasonalProfile,...
          modPnumByVcMonthlyProfile,modPnumByVcSeasonalProfile,...
          uvpPnumMonthlyProfile,uvpPnumSeasonalProfile,uvpPnumMonthlyTargetValues,...
          uvpPnumMonthlyTargetDepths,uvpPnumSeasonalTargetValues,uvpPnumSeasonalTargetDepths] =...
    loadAndFormatModelledParticleNumberDataAndObservations(choiceTypeGridDomain,config,...
        filenameTimeseriesInformation,filenameObsPnumUvp5,fullpathModelRunsDir,filenameSlamsOutput)

% ======================================================================= %
%                                                                         %
% Loads and formats modelled and observational local monthly particle     %
% numebr concentration data in preparation for plotting or model          %
% optimisation.                                                           %
%                                                                         %
% ======================================================================= %

idxSeason = config.seasonIndices;

%% SECTION 1: Load metadata

% Load station-related information from POC flux observation compilation
% (the order of locations is the same as in the UPV5 data)
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
nLocs = length(STATION_NAMES);

% Rearrange model locations to match observational order
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'}; % as in config.gridLats, config.gridLons
[~,reorderModLocIdx] = ismember(STATION_NAMES,currentModLocationOrder); % get reordering indices

%% SECTION 2: Load and process observational data

load(fullfile('.','data','processed',filenameObsPnumUvp5),...
    'uvpPnumByCast','diameterClassEdges','uvpDepths','castMonthlyDistrib')

nEcoTaxaSizeClasses = length(diameterClassEdges);
nObsDepths = length(uvpDepths);

% Calculate the monthly mean
uvpPnumMonthlyProfile = NaN(nEcoTaxaSizeClasses,nObsDepths,12,nLocs);
for iLoc = 1:nLocs   
    for iSc = 1:nEcoTaxaSizeClasses
        for iMonth = 1:12
            nCasts = castMonthlyDistrib(iMonth,iLoc);
            if (nCasts > 0)
                for iDepth = 1:nObsDepths
                    % Number of casts that have reached that depth
                    nCastsInDepth = sum(~isnan(uvpPnumByCast(1:nCasts,iSc,iDepth,iMonth,iLoc)));
                    if (nCastsInDepth > 0)
                        uvpPnumMonthlyProfile(iSc,iDepth,iMonth,iLoc) = mean(uvpPnumByCast(1:nCasts,iSc,iDepth,iMonth,iLoc),'omitnan');
                    end % if nCastsInDepth > 0
                end % iDepth
            end % if nCasts > 0
        end % iMonth
    end % iSc
end % iLoc

% Calculate the seasonal average
uvpPnumSeasonalProfile = NaN(nEcoTaxaSizeClasses,nObsDepths,4,nLocs);
for iSeason = 1:4
    uvpPnumSeasonalProfile(:,:,iSeason,:) = mean(uvpPnumMonthlyProfile(:,:,idxSeason{iSeason},:),3,'omitnan');
end

%% SECTION 3: Load modelled data

load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')

% Assemble modelled particle number array (nSizeClasses x nDepts x 12 months x nLocs)
modPnumByScMonthlyProfile = squeeze(output.particle.avgAttInScMonthly(:,:,1,:,:));  
modPnumByScMonthlyProfile = modPnumByScMonthlyProfile(:,:,:,reorderModLocIdx);

modPnumByVcMonthlyProfile = squeeze(output.particle.avgAttInVcMonthly(:,:,1,:,:));  
modPnumByVcMonthlyProfile = modPnumByVcMonthlyProfile(:,:,:,reorderModLocIdx);

% Calculate seasonal average
modPnumByScSeasonalProfile = NaN(size(modPnumByScMonthlyProfile,1),size(modPnumByScMonthlyProfile,2),4,nLocs);
for iSeason = 1:4
    modPnumByScSeasonalProfile(:,:,iSeason,:) = mean(modPnumByScMonthlyProfile(:,:,idxSeason{iSeason},:),3,'omitnan');
end

modPnumByVcSeasonalProfile = NaN(size(modPnumByVcMonthlyProfile,1),size(modPnumByVcMonthlyProfile,2),4,nLocs);
for iSeason = 1:4
    modPnumByVcSeasonalProfile(:,:,iSeason,:) = mean(modPnumByVcMonthlyProfile(:,:,idxSeason{iSeason},:),3,'omitnan');
end

%% SECTION 4: Extract observed data at modelled depth horizons

[uvpPnumMonthlyTargetValues,uvpPnumMonthlyTargetDepths] = ...
    extractKeyDepthHorizonDataInUvpParticleNumbers(config,uvpPnumMonthlyProfile,uvpDepths);

% Calculate seasonal average
uvpPnumSeasonalTargetValues = NaN(nEcoTaxaSizeClasses,size(uvpPnumMonthlyTargetValues,2),4,nLocs);
uvpPnumSeasonalTargetDepths = NaN(nEcoTaxaSizeClasses,size(uvpPnumMonthlyTargetDepths,2),4,nLocs);
for iSeason = 1:4
    uvpPnumSeasonalTargetValues(:,:,iSeason,:) = mean(uvpPnumMonthlyTargetValues(:,:,idxSeason{iSeason},:),3,'omitnan');
    uvpPnumSeasonalTargetDepths(:,:,iSeason,:) = mean(uvpPnumMonthlyTargetDepths(:,:,idxSeason{iSeason},:),3,'omitnan');
end

end % loadAndFormatModelledParticleNumberDataAndObservations