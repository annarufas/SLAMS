function data = fillNansStrategicallyInSurfaceOceanData(dataOriginal,t,...
    logID,choiceFillingMethod)

% FILLNANSSTRATEGICALLYINSURFACEOCEANDATA Fills temporal gaps in surface 
% ocean time-series data choosing between a standard interp1 interpolation 
% method, or a custom, multi-step interpolation process.
%
% For the custom method, it performs gap-filling in two distinct stages:
% 1. Edge gap: handles gaps at the beginning or end of the time-series 
%    (where valid data exists only on one side) by extrapolating values 
%    from the opposite end using a flip-based interpolation approach.
% 2. Gaps in the middle of the time-series (embedded gaps): interpolates 
%    missing values that are surrounded by valid data on both sides.
%
%   INPUT: 
%       dataOriginal        - data (time vector) containing gaps (NaN values)
%       t                   - time vector indices corresponding to the data (e.g., 1:12 for monthly data)
%       logID               - file to record progress
%       choiceFillingMethod - choose between (0) interp1 and (1) custom (OPTIONAL)
%
%   OUTPUT:
%       data                - data (time vector) with gaps filled in
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

if nargin < 4 || isempty(choiceFillingMethod)
    choiceFillingMethod = 1; % default: custom interpolation
end

% Ensure the input data is a one-dimensional vector
if size(dataOriginal, 2) > 1
    throwError(logID,'ERROR: the data array is not a one dimensonal vector');
end

%% Filling

data = dataOriginal; % make a copy
nanGaps = isnan(data); % identify NaN positions

% Plain interpolation method (interp1)
if choiceFillingMethod == 0  
    data(nanGaps) = interp1(t(~nanGaps),data(~nanGaps),t(nanGaps),'linear','extrap')';

% Custom interpolation method, with two rounds, to handle edge NaNs 
% and embedded NaNs differently    
elseif choiceFillingMethod == 1

    % Round 1: handle edge NaNs (NaNs at the start or end)
    hasEdgeNaNs = isnan(data(1)) || isnan(data(end));
    if hasEdgeNaNs
        data = interpolateEdgeGaps(data,t);
    end  

    % Round 2: handle embedded NaNs (single values or blocks of 
    % various consecutive NaNs that are not at the edges)
    if any(isnan(data))  
        nanGaps = isnan(data); % update
        embeddedNaNs = nanGaps & (1 < (1:length(data)))' & (length(data) > (1:length(data)))';
        data(embeddedNaNs) = interp1(t(~embeddedNaNs),data(~embeddedNaNs),t(embeddedNaNs),'linear','extrap')';
    end

    % Final check: ensure no NaNs are left
    if any(isnan(data))
        throwError(logID,'ERROR: NaNs still left after 2 interpolation rounds');
    end

end % choiceFillingMethod

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

function dataFilled = interpolateEdgeGaps(dataOriginal,t)
 
    % This function uses a flip-based interpolation method
    
    % Make a copy
    dataFilled = dataOriginal;

    % Count NaNs from the beginning and the end
    nStartNans = find(~isnan(dataOriginal), 1) - 1; 
    nEndNans = length(dataOriginal) - find(~isnan(dataOriginal), 1, 'last'); 
    
    % Determine the number of edges that need to be filled in (either "1"
    % or "2")
    nNanEdges = (nStartNans > 0) + (nEndNans > 0);
    
    % Identify the first and last valid values
    idxFirstValid = find(dataOriginal >= 0,1,'first'); 
    idxLastValid = find(dataOriginal >= 0,1,'last'); 
    
    % Prepare interpolation vector
    lengthGap = nStartNans + nEndNans;
    wrapGapData = NaN(lengthGap+2, 1);
    wrapGapData(1) = dataOriginal(idxFirstValid);
    wrapGapData(end) = dataOriginal(idxLastValid);
    matchValidIdxs = find(~isnan(wrapGapData));
    matchNanIdxs = find(isnan(wrapGapData));
    
    % Vector filled in
    gapFilled = interp1(matchValidIdxs,wrapGapData(matchValidIdxs),...
        matchNanIdxs,'linear','extrap');
    
    % Handle the filling based on the number of edges
    if nNanEdges == 1 % either at the beginning or end, i.e., #######-------- (or) -------#######
     
        if idxFirstValid > t(1) % fill the start of the time series 
            dataFilled(1:idxFirstValid-1) = flip(gapFilled); 
        elseif idxFirstValid == t(1) % fill the end of the time series 
            dataFilled(idxLastValid+1:end) = flip(gapFilled); 
        end   
   
    elseif nNanEdges == 2 % both at the beginning and end, i.e.,  ------#####----- 
        
        % Fill both ends of the time series 
        nValsFirstSection = idxFirstValid-1;
        dataFilled(1:nValsFirstSection) = flip(gapFilled(1:nValsFirstSection));
        dataFilled(idxLastValid+1:end) = flip(gapFilled(nValsFirstSection+1:end));

    end
             
end % interpolateEdgeGaps
   
end % fillNansStrategicallyInSurfaceOceanData