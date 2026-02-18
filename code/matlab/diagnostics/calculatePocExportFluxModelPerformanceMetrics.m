function metricScores = calculatePocExportFluxModelPerformanceMetrics(...
    obsFluxDataMonthly,bicepExFluxMonthly,modFluxMonthlyDh,...
    filenameTimeseriesInformation,pathToOutputFile)

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
nLocs = length(STATION_NAMES);

% Extract data at POC export flux where necessary
modPocExport = squeeze(modFluxMonthlyDh(1,:,:,1,1));
insituPocExport = squeeze(obsFluxDataMonthly{1}(1,:,:)); 
satellPocExport = permute(bicepExFluxMonthly, [2, 3, 1]); % from nLocs x 12 x nModels --> 12 x nModels x nLocs

% Rearrange locations to match desired order
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
currentLocationOrder = STATION_NAMES;
[~,reorderLocIdx] = ismember(desiredLocationOrder,currentLocationOrder); % get reordering indices
nLocs = length(desiredLocationOrder);

modPocExport = modPocExport(:,reorderLocIdx);
insituPocExport = insituPocExport(:,reorderLocIdx);
satellPocExport = satellPocExport(:,:,reorderLocIdx);

% Organise data into cell arrays
obsAll = {insituPocExport, squeeze(satellPocExport(:,1,:)), squeeze(satellPocExport(:,2,:)), squeeze(satellPocExport(:,3,:))};
modAll = {modPocExport, modPocExport, modPocExport, modPocExport};

nVars = length(obsAll);
labelVars = {'In situ','BICEP-Dunne','BICEP-Henson','BICEP-Li'}; 

% Metric calculation
nMetrics = 6; % number of statistics to compute
metricScores = NaN(nLocs,nVars,nMetrics);

for iVar = 1:nVars
    obsVarData = obsAll{iVar};
    modVarData = modAll{iVar};

    for iLoc = 1:nLocs
        localObs = obsVarData(:,iLoc); 
        localMod = modVarData(:,iLoc);

        % Filter out NaN observations and corresponding modelled data
        validIdx = ~isnan(localObs) & ~isnan(localMod);
        localObsValid = localObs(validIdx);
        localModValid = localMod(validIdx);

        % Compute matchup statistics
        [metricScores(iLoc,iVar,:),labelMetrics] = calculatePerformanceMetrics(...
            localObsValid,localModValid);
    end
end

% Generate LaTeX table for regression statistics
generateMetricsLatexOutput(metricScores,desiredLocationOrder,labelVars,...
    labelMetrics,pathToOutputFile);

end % calculatePocExportFluxModelPerformanceMetrics