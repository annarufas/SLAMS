function data = fillNansStrategicallyInDepthDependentOceanData(dataOriginal,logID)

% FILLNANSSTRATEGICALLYINDEPTHDEPENDENTOCEANDATA Fills NaN gaps over depth 
% and time for a specific location in depth-dependent ocean time-series 
% data (e.g., temperature, nutrients, oxygen) using strategic interpolation.
%
%   INPUT: 
%       dataOriginal - data array (depth x time) containing gaps (NaN values)
%       logID        - file to record progress
%
%   OUTPUT:
%       data         - data array (depth x time) with gaps filled in
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 1 Feb 2025 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Assumption check

% Create mask for valid values in the first time step. We assume that the 
% rest of time steps follow same depth pattern as season does not affect 
% the vertical extent of data, unlike surface ocean data)
validPointsFirstTimeStep = ~isnan(dataOriginal(:,1));

% Check above assumption and ensure NaN patterns are consistent across time 
% steps
nanPattern = isnan(dataOriginal);
consistentPattern = all(nanPattern == nanPattern(:, 1), 2); % compare every column with the first column
if ~all(consistentPattern)
    throwError(logID, 'ERROR: the assumption that the NaN pattern remains constant across time steps is incorrect for this dataset.\n');
end

%% Filling NaNs

data = dataOriginal; % create a copy
[~, nTimeSteps] = size(dataOriginal);

if any(validPointsFirstTimeStep)
    for iTimeStep = 1:nTimeSteps
        depthSlice = dataOriginal(:,iTimeStep);

        % Detect the first and last valid values
        firstValidIdx = find(~isnan(depthSlice), 1, 'first');
        lastValidIdx = find(~isnan(depthSlice), 1, 'last');

        % Edge NaN indices
        edgeNaNsStart = 1:firstValidIdx - 1;
        edgeNaNsEnd = lastValidIdx + 1:length(depthSlice);

        % Detect embedded NaNs (ignore the edge ones)
        innerData = depthSlice(firstValidIdx:lastValidIdx);
        embeddedNaNs = find(isnan(innerData)) + firstValidIdx - 1;

        if ~isempty(embeddedNaNs)
            z = 1:length(depthSlice);
            validMask = ~isnan(depthSlice);
            if sum(~isnan(depthSlice))
                depthSlice(embeddedNaNs) = interp1(z(validMask),...
                    depthSlice(validMask),z(embeddedNaNs),'linear','extrap')';
            end
        end 
        if ~isempty(edgeNaNsEnd) 
            depthSlice(edgeNaNsEnd) = depthSlice(lastValidIdx); % forward-fill edge NaNs at the end
        end  
        if ~isempty(edgeNaNsStart)
            depthSlice(edgeNaNsStart) = depthSlice(firstValidIdx); % backward-fill edge NaNs at the start
        end 

        data(:,iTimeStep) = depthSlice;

    end
end
    
end % fillNansStrategicallyInDepthDependentOceanData
