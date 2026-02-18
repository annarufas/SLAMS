function metricScores = calculateOverallModelPerformanceMetrics(config,...
    obsFluxDataMonthly,modFluxMonthlyDh,uvpPnumMonthlyTargetValues,...
    modPnumMonthlyProfile,uvpDepths,filenameTimeseriesInformation,pathToOutputFile)

% CALCULATEOVERALLMODELPERFORMANCEMETRICS Computes performance metrics 
% between modelled outputs and corresponding observations and generates a 
% .tex table for outputs.

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
nLocs = length(STATION_NAMES);

% Rearrange flux model outputs and observations: (12*nMatchedDepths x nTracers x nLocs)
nRows = size(obsFluxDataMonthly{1},1) * size(obsFluxDataMonthly{1},2);
arrayFluxObs = NaN(nRows,3,nLocs);
arrayFluxMod = NaN(nRows,3,nLocs);
for iTracer = 1:3
    for iLoc = 1:nLocs
        tracerObs = squeeze(obsFluxDataMonthly{iTracer}(:,:,iLoc)); 
        tracerMod = squeeze(modFluxMonthlyDh(:,:,iLoc,1,iTracer));
        arrayFluxObs(:,iTracer,iLoc) = reshape(tracerObs,[],1);
        arrayFluxMod(:,iTracer,iLoc) = reshape(tracerMod,[],1); 
    end
end

% Before rearranging pnum, sum over particle size classes and get equal
% number of depths
upvTotPnum = squeeze(sum(uvpPnumMonthlyTargetValues,1,'omitnan'));
modTotPnum = squeeze(sum(modPnumMonthlyProfile,1,'omitnan'));

% Cap modelledDepths where uvpDepths stop
modelledPnumDepths = config.availImagSysDeployDepths;
[~,idxLastPnumDepth] = min(abs(modelledPnumDepths - uvpDepths(end)));
modelledPnumDepths = modelledPnumDepths(1:idxLastPnumDepth);
modTotPnum = modTotPnum(1:length(modelledPnumDepths),:,:);

% Rearrange particle num. concentration model outputs and observations (12*nMatchedDepths x nLocs)
nRows = size(upvTotPnum,1) * size(upvTotPnum,2);
arrayPnumObs = NaN(nRows,nLocs);
arrayPnumMod = NaN(nRows,nLocs);
for iLoc = 1:nLocs
    pnumObs = squeeze(upvTotPnum(:,:,iLoc)); 
    pnumMod = squeeze(modTotPnum(:,:,iLoc));
    arrayPnumObs(:,iLoc) = reshape(pnumObs,[],1);
    arrayPnumMod(:,iLoc) = reshape(pnumMod,[],1); 
end

% Rearrange locations to match desired order
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
currentLocationOrder = STATION_NAMES;
[~,reorderLocIdx] = ismember(desiredLocationOrder,currentLocationOrder); % get reordering indices
nLocs = length(desiredLocationOrder);

arrayFluxObs = arrayFluxObs(:,:,reorderLocIdx);
arrayFluxMod = arrayFluxMod(:,:,reorderLocIdx);
arrayPnumObs = arrayPnumObs(:,reorderLocIdx);
arrayPnumMod = arrayPnumMod(:,reorderLocIdx);

% Organise data into cell arrays
obsAll = {squeeze(arrayFluxObs(:,1,:)), squeeze(arrayFluxObs(:,2,:)), squeeze(arrayFluxObs(:,3,:)), arrayPnumObs};
modAll = {squeeze(arrayFluxMod(:,1,:)), squeeze(arrayFluxMod(:,2,:)), squeeze(arrayFluxMod(:,3,:)), arrayPnumMod};

nVars = length(obsAll);
labelVars = {'POC flux','PIC flux','bSi flux','Particle num. concentration'}; 

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

end % calculateOverallModelPerformanceMetrics