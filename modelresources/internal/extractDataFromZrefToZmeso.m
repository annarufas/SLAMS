function [selectedDepths,selectedFluxes,selectedErrs] = extractDataFromZrefToZmeso(...
    depthVector,fluxVector,fluxErrVector,choiceZref,depthBoundaries,maxZref)

% EXTRACTDATAFROMZREFTOZMESO Extracts POC flux values and associated depths 
% below the reference depth (zref) up to the base of the emsopelagic zone 
% (zmeso).
%
%   INPUT: 
%       depthVector     - vector of depth values
%       fluxVector      - vector of POC flux values
%       fluxErrVector   - vector of POC flux errors
%       choiceZref      - choice of reference depth (zeu, 100 m or inflexion point)
%       depthBoundaries - local euphotic layer depth, mesopelagic layer
%                         depth and enar seafloor depth boundaries
%       maxZref         - maximum allowable reference depth
%
%   OUTPUT:
%       selectedDepths - vector of selected depths
%       selectedFluxes - vector of selected POC flux values
%       selectedErrs   - vector of selected POC flux errors
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 25 Nov 2024  
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------
                
%% Get rid of NaNs

validIdx = ~isnan(depthVector) & ~isnan(fluxVector); 
validDepths = depthVector(validIdx);
validFluxes = fluxVector(validIdx);
validErrs = fluxErrVector(validIdx);

%% Initialise outputs

selectedDepths = [];
selectedFluxes = [];
selectedErrs = [];

%% Determine zref (based on choiceZref) and zmeso

% Reference depth
switch choiceZref
    case 1 % Fixed 100 m depth
        [~,idxZref] = min(abs(validDepths - 100)); 
        
    case 2 % Local euphotic depth (zeu)
        idxZref = [];
        
        % Check if depths are between the boundaries
        inRange = validDepths >= depthBoundaries(1,1) & validDepths <= depthBoundaries(2,1);

        % Localise the closest depth to zeu
        if sum(inRange) > 0
            locZeu = (depthBoundaries(2,1) + depthBoundaries(1,1))/2;
            [~,idxZref] = min(abs(validDepths - locZeu)); 
        end
        
    case 3 % Depth at inflection point 
        % Find the depth where flux first decreases (simple check for inflection)
        idxZref = [];
        for i = 1:length(validFluxes)-1
            if all(validFluxes(i) > validFluxes(i+1:end))
                idxZref = i;
                break;
            end
        end

end

% We have chosen to fit b and z* up to zmeso, so that is the depth we are
% goig to extract here
idxZmeso = [];
% Check if depths are between the boundaries
inRange = validDepths >= depthBoundaries(1,2) & validDepths <= depthBoundaries(2,2);
% Localise the closest depth to zeu
if sum(inRange) > 0
    locZmeso = (depthBoundaries(2,2) + depthBoundaries(1,2))/2;
    [~,idxZmeso] = min(abs(validDepths - locZmeso)); 
end

%% Extract

if ~isempty(idxZref) && ~isempty(idxZmeso)
    depthCondition = validDepths >= validDepths(idxZref) & validDepths <= validDepths(idxZmeso);
    selectedDepths = validDepths(depthCondition);
    
    % Ensure that selected depths do not exceed maxZref
    if ~isempty(selectedDepths) && selectedDepths(1) <= maxZref
        selectedFluxes = validFluxes(depthCondition); 
        selectedErrs = validErrs(depthCondition); 
    else
        selectedDepths = [];
        selectedFluxes = [];
        selectedErrs = [];
    end
end

end % extractDataFromZrefToZmeso