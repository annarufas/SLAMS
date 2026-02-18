function dataFilled = manageGapsInOceanData(dataOriginal,z,t,hasDepth,F,...
    currentLat,currentLon,lonIncrement,latIncrement,lastLat,maxIterations,logID)
        
% MANAGEGAPSINOCEANDATA Fills gaps over depth and time for a
% specific location in depth-dependent ocean time-series data (e.g., 
% temperature, nutrients, oxygen) using strategic interpolation and 
% neighborhood averaging techniques. 
%
% This function manages missing data by:
% 1. Filling isolated NaNs through interpolation.
% 2. Expanding the spatial search if the entire time-depth series at a location is NaN. 
% 3. Attempting final interpolation if partial gaps remain after neighborhood search.
% 
%   INPUT: 
%       dataOriginal  - data array (depth x time or just time) containing NaN values at a specific location lat/lon
%       z             - vector of depth levels corresponding to dataOriginal rows
%       t             - time vector indices corresponding to the data (e.g., 1:12 for monthly data)
%       hasDepth      - boolean indicating whether the data has depth dimension
%       F             - gridded data interpolant used for neighborhood averaging
%       currentLat    - current latitude of the data point
%       currentLon    - current longitude of the data point
%       lonIncrement  - longitude increment for expanding the search area
%       latIncrement  - latitude increment for expanding the search area
%       lastLat       - boundary latitude to constrain the search area
%       maxIterations - maximum number of iterations for spatial search
%       logID         - file to record progress
%
%   OUTPUT:
%       dataFilled    - data array (depth x time or just time) with gaps filled in
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

%% Strategic filling

dataFilled = dataOriginal; % make a copy

% Perform strategic filling only if:
% 1) There are NaN values present in the data.
% 2) Not all values are NaN (i.e., it is not a land cell with entirely missing data).
if any(isnan(dataOriginal(:))) && ~all(isnan(dataFilled(:)))
    % 2D depth-dependent data (no minimum valid values required)
    if hasDepth
    	%disp('has depth')
        dataFilled = fillNansStrategicallyInDepthDependentOceanData(dataOriginal,logID);
    % 1D surface ocean data (requires at least 2 valid values due to interpolation method)
    elseif ~hasDepth && sum(~isnan(dataOriginal(:))) > 2 
        dataFilled = fillNansStrategicallyInSurfaceOceanData(dataOriginal,t,logID);
    end
end

%% Neighborhood mean search for completely missing data

% This section addresses cases where all values (depth and time) at a specific 
% coordinate point are NaN. These scenarios occur when initial strategic 
% filling failed to yield any valid data.

if any(isnan(dataFilled(:))) && ~all(isnan(dataFilled(:)))
    
    % Error handling: Sparse NaNs remaining unexpectedly
    % This situation should not occur—data should either have no NaNs or 
    % be completely filled with NaNs. Remaining isolated NaNs indicate an error.
    fprintf(logID,'WARNING: sparse NaNs in dataset, only one data point. Addressing...\n');
    
    % Fill by repeating the only value present as we cannot leave NaN
    % positions
    nanGaps = isnan(dataFilled); % identify NaN positions
    dataFilled(nanGaps) = dataFilled(~isnan(dataFilled)); % fill with only valid value

elseif all(isnan(dataFilled(:))) 
	%disp('all is nan')
    % Determine whether the data has depth levels or is surface-only
    if hasDepth
        depthLevels = z;
        %disp('depth levels')
        %disp(numel(depthLevels))
    else
        depthLevels = [];  % surface data has no depth levels
    end

    % Search for valid neighborhood data and compute spatial mean
    for iDepth = 1:max(1,numel(depthLevels))  
        currentDepth = depthLevels; 
        if ~isempty(depthLevels)
            currentDepth = depthLevels(iDepth);
        end
        %disp('iDepth')
        %disp(iDepth)
        %disp(currentDepth)
        areaAverage = findNeighbourhoodMean(F,currentLon,currentLat,lonIncrement,...
            latIncrement,lastLat,maxIterations,currentDepth,t,logID);
        if ~isempty(depthLevels)
            dataFilled(iDepth,:) = areaAverage;
        else
            dataFilled = areaAverage;
        end
    end 

    % Attempt another strategic fill if NaNs remain
    if any(isnan(dataFilled(:))) && ~all(isnan(dataFilled(:)))   
        if hasDepth
            fprintf(logID,'Some NaNs left after neighborhood search, attempting to fill them.\n');
        	dataFilled = fillNansStrategicallyInDepthDependentOceanData(dataFilled,logID);
            if any(isnan(dataFilled(:)))
                throwError(logID,'ERROR: filling was not possible.\n');
            end
        elseif ~hasDepth && sum(~isnan(dataFilled(:))) > 2
            fprintf(logID,'Some NaNs left after neighborhood search, attempting to fill them.\n');
            dataFilled = fillNansStrategicallyInSurfaceOceanData(dataFilled,t,logID);
            if any(isnan(dataFilled(:)))
                throwError(logID,'ERROR: filling was not possible.\n');
            end
        end
    end    
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function areaAverage = findNeighbourhoodMean(F,startLon,startLat,...
    lonIncrement,latIncrement,lastLat,maxIterations,z,t,logID)

    % This function searches iteratively for a neighborhood containing 
    % valid data in a grid dataset. If valid values are found, it computes 
    % the average over time for the given spatial region. Otherwise, it 
    % expands the search area by increasing the latitude and longitude 
    % boundaries until valid data is found or the maximum number of 
    % iterations is reached. If no valid neighborhood is found after all 
    % iterations, the function returns NaN for the average.

    stepSizeLon = lonIncrement; % longitude increment for neighborhood expansion
    stepSizeLat = latIncrement; % latitude increment for neighborhood expansion
    
    areaAverage = NaN(length(t),1); % initialise output array
    iIteration = 1; % initialise search iteration counter
    
    % It should be possible to find values within the maxIterations that we
    % have set as those are coastal cases in the viccinity of open ocean
    % data
    while iIteration <= maxIterations

        % Calculate boundary coordinates for current search region
        latUp = min(startLat + latIncrement, lastLat(end)); % esure latitude stays within valid range
        latDown = max(startLat - latIncrement, -1*lastLat(1)); % same as above

        lonRight = startLon + lonIncrement; 
        lonLeft = startLon - lonIncrement;
        lonRight = mod(lonRight + 180, 360) - 180; % ensure longitude stays within [-180, 180]
        lonLeft = mod(lonLeft + 180, 360) - 180; % same as above
        
        % Generate latitude and longitude search ranges
        latRange = latDown:stepSizeLat:latUp;
        if lonLeft < lonRight
            % Simple case: longitude range does not cross the 180º meridian
            lonRange = lonLeft:stepSizeLon:lonRight;
        else
            % Complex case: longitude range crosses the 180º meridian
            % Compute the edge point past 180º to seamlessly continue the range
            edge = -180 + abs(180 - (lonLeft + stepSizeLon)); % to resolve problems when crossing the 180 parallel
            lonRange = [lonLeft:stepSizeLon:180, edge:stepSizeLon:lonRight];
        end
        
        % Log search ranges for debugging
        if isempty(z) || z < 7 % only print info for first depth
            fprintf(logID, 'In iteration %d, lat range: [%.2f to %.2f], lon range: [%.2f to %.2f]\n', ...
                iIteration, latDown, latUp, lonLeft, lonRight);
        end
        
        % Extract neighborhood data by generating a 4D grid (lat, lon, z, time)
        if isempty(z)
            [qX, qY, qT] = ndgrid(latRange,lonRange,t');
            interpolatedData = squeeze(F(qX, qY, qT)); % 3D array (nLat x nLon x time)
        else
            [qX, qY, qZ, qT] = ndgrid(latRange,lonRange,z',t');
            interpolatedData = squeeze(F(qX, qY, qZ, qT)); % 3D array (nLat x nLon x time)
        end
        %disp('size interpolated data')
        %disp(size(interpolatedData))

        % Compute the mean of valid data points, ignoring NaN values
        if sum(~isnan(interpolatedData(:))) > 0
        %disp('numel t')
        %disp(numel(t))
        %disp('t')
        %disp(t)
            for iIdx = 1:numel(t)
                timeSliceForTheArea = interpolatedData(:,:,iIdx);
                areaAverage(iIdx) = mean(timeSliceForTheArea(:), 'omitnan');
            end
            if isempty(z) || z < 7 % only print info for first few depths
                fprintf(logID, 'Mean found\n');
            end
            return;
        else
            if isempty(z) || z < 7 % only print info for first few depths
                fprintf(logID, 'No valid neighbours found after a full search at depth %d\n',z);
            end
        end
        
        % Update increments for next search iteration
        lonIncrement = lonIncrement + stepSizeLon;
        latIncrement = latIncrement + stepSizeLat;
        iIteration = iIteration + 1;
        
    end

end % findNeighbourhoodMean

% *************************************************************************
    
end % manageGapsInOceanData
