function data = manageGapsInOceanDataIceMasked(dataOriginal,mask,...
    data_lat,data_lon,mask_lat,mask_lon,t,logID,maskThreshold,choiceFillingMethod)

% MANAGEGAPSINOCEANDATAICEMASKED Fills temporal gaps in surface ocean
% time-series data affected by ice coverage (e.g., chla, PAR0, NPP, dust 
% flux) choosing between a standard interp1 interpolation method, or a 
% custom, multi-step interpolation process. After interpolation, it applies 
% an ice-based mask, which will set to NaN grid cells that do not meet a 
% mask threshold criteria. Notice that this fucntion is for surface ocean
% data only, without depth-dependency.
%
% For the custom method, it performs gap-filling in two distinct stages:
% 1. Edge gap: handles gaps at the beginning or end of the time-series 
%    (where valid data exists only on one side) by extrapolating values 
%    from the opposite end using a flip-based interpolation approach.
% 2. Gaps in the middle of the time-series (embedded gaps): interpolates 
%    missing values that are surrounded by valid data on both sides.
%
%   INPUT: 
%       dataOriginal        - data array (time) containing gaps (NaN values)
%       mask                - mask array
%       data_lat            - latitude vector for data 
%       data_lon            - longitude vector for data
%       mask_lat            - latitude vector for the mask 
%       mask_lon            - longitude vector for the mask
%       t                   - time vector indices corresponding to the data (e.g., 1:12 for monthly data)
%       logID               - file to record progress
%       maskThreshold       - threshold above which ocean colour products become NaN (OPTIONAL)
%       choiceFillingMethod - choose between (0) interp1 and (1) custom (OPTIONAL)
%
%   OUTPUT:
%       data                - data array (time) with gaps filled in and mask applied
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

%% Input validation

% Check if data and mask dimensions match
if ndims(dataOriginal) ~= ndims(mask)
    throwError(logID,'ERROR: input data array and mask must have the same dimensions');
end

% Check time vector length based on dimensions
if ndims(dataOriginal) == 3 % a lat x lon x time array
    if length(t) ~= size(dataOriginal, 3)
        throwError(logID,'ERROR: time vector length must match the 3rd dimension of the data array');
    end
else
    throwError(logID,'ERROR: the data array does not have dimensions lat x lon x time');
end

% Handle optional arguments
if nargin < 10 || isempty(choiceFillingMethod)
    choiceFillingMethod = 1; % default: custom interpolation
end
if nargin < 9 || isempty(maskThreshold)
    maskThreshold = []; % default: use with no threshold
end

%% Regrid the mask to same grid used by data

if ndims(dataOriginal) == 3 && ~isequal(size(dataOriginal), size(mask))
    mask = regridVariable(data_lat,data_lon,[],t,mask,mask_lat,mask_lon,[],t);  
end

%% Filling and masking

data = dataOriginal; % copy the dataset
[nLat,nLon,~] = size(dataOriginal);

for iRow = 1:nLat % latitudes
    for iCol = 1:nLon % longitudes

        % Extract the time series at the current grid point
        localData = squeeze(dataOriginal(iRow,iCol,:)); 
        localMask = squeeze(mask(iRow,iCol,:));
        
        % Perform interpolation only if:
        % 1) There are NaN values present.
        % 2) Not all values are NaN (i.e., this is not a permanent land cell).
        if any(isnan(localData)) && sum(~isnan(localData)) > 2
            data(iRow,iCol,:) = fillNansStrategicallyInSurfaceOceanData(localData,...
                t,logID,choiceFillingMethod);
        end
        
        % Apply mask
        data(iRow,iCol,:) = applySurfaceOceanMask(squeeze(data(iRow,iCol,:)),...
            localData,localMask,maskThreshold);
    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function timeseriesFilled = applySurfaceOceanMask(timeseriesFilled,timeseriesOriginal,...
    timeseriesMask,maskThreshold)
    
    % If no mask threshold is provided, assume binary mask (0: masked, 1: valid)
    if isempty(maskThreshold)
        % Set to NaN cells where mask = 0
        maskInvalidPoints = (timeseriesMask == 0);
        timeseriesFilled(maskInvalidPoints) = NaN;

    % Use the mask threshold value provided    
    else  
        % Set NaN only where:
        % - the mask value exceeds the threshold, AND
        % - the corresponding point in the original time series is already NaN
        % This avoids unnecessarily discarding valid underlying data.
        maskInvalidPoints = (timeseriesMask > maskThreshold) & isnan(timeseriesOriginal);
        timeseriesFilled(maskInvalidPoints) = NaN;
    end

end % applySurfaceOceanMask

end % manageGapsInOceanDataIceMasked
