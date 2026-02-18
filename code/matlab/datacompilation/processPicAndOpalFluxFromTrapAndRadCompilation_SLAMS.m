
% ======================================================================= %
%                                                                         %
% This script reads in the PIC and bSi flux compilation, calculates       %
% monthly and annual averages and propagates error accordingly. It is     %
% based on the script to analyse POC flux data that I created for the     %
% 2025 GRL publication, processPocFluxFromTrapAndRadCompilation.m, which  %
% only needed a few code additions to consider two variables instead of   %
% one. Output units for PIC flux are "mmol C m-2 d-1", and for bSi flux   %
% "mmol Si m-2 d-1".                                                      %
%                                                                         %
% The script has 8 sections:                                              %
%   Section 1 - Presets.                                                  %
%   Section 2 - Load the dataset and manipulate the data array.           %
%   Section 3 - Bin data monthly by depth horizon and propagate error.    %
%   Section 4 - Bin data monthly by unique depth and propagate error.     %
%   Section 5 - Bin data annually by depth horizon and propagate error.   %
%   Section 6 - Bin data annually by unique depth and propagate error.    %
%   Section 7 - Calculate the number of data points based on various      % 
%               criteria.                                                 %
%   Section 8 - Save the data.                                            %              
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 8 May 2025                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/raw/observations/'));
addpath(genpath('./data/processed/'));
addpath(genpath('./code/'));
addpath(genpath('./modelresources/external/'));
addpath(genpath('./modelresources/internal/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Filename declarations 
filenameInputFluxCompilation       = 'dataset_s0_pic_bsi_flux_trap_and_radionuclide_compilation.xlsx';
filenameInputZeu                   = 'zeu_calculated_chlaoccci_mldifremer_pointonepercentpar0.mat';
filenameInputTimeseriesInformation = 'timeseries_station_information_slams.mat';
filenameOutputFluxCompilation      = 'picandbsiflux_compilation_slams.mat';

% Define parameters
RAND_ERR_FRAC_PIC = 0.30; % 30% (Buesseler et al. 2000, Buesseler et al. 2007, Stanley et al. 2004)
RAND_ERR_FRAC_SI = 0.30;
SYS_ERR_FRAC_PIC = 0.05;
SYS_ERR_FRAC_SI = 0.075; 
OC_ERR_FRAC = 0.10; % ocean colour error fraction (10%), based on McKinna et al. 2019
MOLAR_MASS_CARBON = 12.011; % g mol-1
MOLAR_MASS_OPAL = 67.3; % g mol-1

% Load in stations information (already created by the script
% processPocFluxFromTrapAndRadCompilation.m)
load(fullfile('.','data','interim',filenameInputTimeseriesInformation),...
    'LOC_LATS','LOC_LONS','LOC_DEPTH_HORIZONS','STATION_NAMES','STATION_TAGS',...
    'MAX_NUM_VALUES_PER_MONTH','MAX_NUM_DEPTHS_PER_PROFILE','qZeuMonthly','MAX_ZEU')

tagDepthHorizons = {'zeu','mesoupp','mesolow','zmeso'};
NUM_TARGET_DEPTHS = length(tagDepthHorizons);
NUM_LOCS = length(STATION_NAMES);

% Duplicate and adapt depth horizons along 5th dimension
LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT = cat(5, LOC_DEPTH_HORIZONS, LOC_DEPTH_HORIZONS);

% HOT/ALOHA site: adjust POC trap target depth (PIC only)      
LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(:,5,1,3,1) = 490;   
LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(:,5,2,3,1) = 510;

% Replicate to fill later
LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD THE DATASET AND MANIPULATE THE DATA ARRAY
% -------------------------------------------------------------------------

% Load the excel spreadhseet with flux data
opts = detectImportOptions(fullfile('.','data','raw',filenameInputFluxCompilation));
opts = setvartype(opts,{'PIC_mmol_m2_d','PIC_mg_m2_d','randerr_PIC_mmol_m2_d','randerr_PIC_mg_m2_d','Si_mmol_m2_d','bSi_mg_m2_d','randerr_Si_mmol_m2_d','randerr_bSi_mg_m2_d'},'double');
opts = setvartype(opts,{'deploymentDate','midDate','recoveryDate'},'datetime');
opts = setvartype(opts,{'method','source','originalUnits','comments'},'string');
D = readtable(filenameInputFluxCompilation,opts);

% Add a 'month' and a 'year' column
D.month = month(D.midDate);
D.year = year(D.midDate);

% Transform into categorical the variable 'month' and 'tag'
D.tag = categorical(D.tag);
D.month = categorical(D.month);

% Add a column to indicate the depth horizon of the data 
D.depthHorizon = cell(height(D),1);

% Make sure the timetable is sorted by 'tag' - ALPHABETICAL ORDER
D = sortrows(D,'tag');

% Add the sys error column
D.syserr_PIC_mmol_m2_d = SYS_ERR_FRAC_PIC.*D.PIC_mmol_m2_d;
D.syserr_Si_mmol_m2_d = SYS_ERR_FRAC_SI.*D.Si_mmol_m2_d;

% Manipulate the random error
% Calculate random error in mmol/m^2/d and update if necessary
D.randerr_PIC_mmol_m2_d(isnan(D.randerr_PIC_mmol_m2_d)) = RAND_ERR_FRAC_PIC .* D.PIC_mmol_m2_d(isnan(D.randerr_PIC_mmol_m2_d));
D.randerr_Si_mmol_m2_d(isnan(D.randerr_Si_mmol_m2_d)) = RAND_ERR_FRAC_SI .* D.Si_mmol_m2_d(isnan(D.randerr_Si_mmol_m2_d));

% Update random error in mg/m^2/d
D.randerr_PIC_mg_m2_d = MOLAR_MASS_CARBON .* D.randerr_PIC_mmol_m2_d;
D.randerr_bSi_mg_m2_d = MOLAR_MASS_OPAL .* D.randerr_Si_mmol_m2_d;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - BIN DATA MONTHLY BY DEPTH HORIZON AND PROPAGATE ERROR
% -------------------------------------------------------------------------

% The variables that we want to extract
tagPicValues = 'PIC_mmol_m2_d';
tagPicSysError = 'syserr_PIC_mmol_m2_d';
tagPicRandError = 'randerr_PIC_mmol_m2_d';
tagSiValues = 'Si_mmol_m2_d';
tagSiSysError = 'syserr_Si_mmol_m2_d';
tagSiRandError = 'randerr_Si_mmol_m2_d';
NUM_VARS = 2;

% Define output arrays
picbsiFluxRawProfileValues_cell   = cell(12,NUM_LOCS,NUM_VARS); 
picbsiFluxRawProfileErrRand_cell  = cell(12,NUM_LOCS,NUM_VARS); 
picbsiFluxRawProfileErrSys_cell   = cell(12,NUM_LOCS,NUM_VARS);
picbsiFluxRawProfileDepths_cell   = cell(12,NUM_LOCS,NUM_VARS); 
picbsiFluxRawProfileDataType_cell = cell(12,NUM_LOCS,NUM_VARS); 

picbsiFluxRawDhValues_cell   = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS,NUM_VARS); 
picbsiFluxRawDhDepths_cell   = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS,NUM_VARS); 
picbsiFluxRawDhDataType_cell = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS,NUM_VARS); % sediment trap vs radionuclides
picbsiFluxRawDhTag_cell      = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS); % zeu, mesoupp, mesolow or zmeso

picbsiFluxMonthlyDhAvg     = NaN(NUM_TARGET_DEPTHS,12,NUM_LOCS,NUM_VARS); % mean
picbsiFluxMonthlyDhN       = zeros(NUM_TARGET_DEPTHS,12,NUM_LOCS,NUM_VARS); % number of values
picbsiFluxMonthlyDhErrTot  = NaN(NUM_TARGET_DEPTHS,12,NUM_LOCS,NUM_VARS); % total error

for iVar = 1:NUM_VARS
    
    if (iVar == 1)
        tagValues = tagPicValues;
        tagSysError = tagPicSysError;
        tagRandError = tagPicRandError;
    elseif (iVar == 2)
        tagValues = tagSiValues;
        tagSysError = tagSiSysError;
        tagRandError = tagSiRandError;
    end
    
    for iLoc = 1:NUM_LOCS
        currStatData = D(D.tag == STATION_NAMES{iLoc},:);

        % Locate the variables that are relevant
        colIdsAll = zeros(6,1);
        colIdsAll(1) = find(strcmpi(currStatData.Properties.VariableNames,tagValues));
        colIdsAll(2) = find(strcmpi(currStatData.Properties.VariableNames,tagSysError));
        colIdsAll(3) = find(strcmpi(currStatData.Properties.VariableNames,tagRandError));
        colIdsAll(4) = find(strcmpi(currStatData.Properties.VariableNames,'month'));
        colIdsAll(5) = find(strcmpi(currStatData.Properties.VariableNames,'depth'));
        colIdsAll(6) = find(strcmpi(currStatData.Properties.VariableNames,'method'));

        % Extract these variables
        currStatData = currStatData(:,colIdsAll');
        
        % Identify rows with any NaN values
        rowsWithNaN = ismissing(currStatData(:,1)); %any(ismissing(currStatData), 2);

        % Remove rows with NaN values
        currStatDataValid = currStatData(~rowsWithNaN, :);

        % (1) Pull all the data at that station and save it to output arrays
        for iMonth = 1:12
            currMonthData = currStatDataValid(currStatDataValid.month == num2str(iMonth),:);
            if (~isempty(currMonthData))
                picbsiFluxRawProfileValues_cell{iMonth,iLoc,iVar} = num2str(currMonthData{:,1}'); 
                picbsiFluxRawProfileErrSys_cell{iMonth,iLoc,iVar} = num2str(currMonthData{:,2}');
                picbsiFluxRawProfileErrRand_cell{iMonth,iLoc,iVar} = num2str(currMonthData{:,3}'); 
                picbsiFluxRawProfileDepths_cell{iMonth,iLoc,iVar} = num2str(currMonthData{:,5}');
                picbsiFluxRawProfileDataType_cell{iMonth,iLoc,iVar} = strjoin(currMonthData{:,6}'); 
            end  
        end

        % (2) Pull data by month
        for iMonth = 1:12
            monthCondition = currStatDataValid.month == num2str(iMonth);
            currStatMonthData = currStatDataValid(monthCondition,:);
      
            % (3) Pull data by the depth horizon of interest
            for iDh = 1:NUM_TARGET_DEPTHS

                % Extract original depth bounds for the layer
                zLower = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar);
                zUpper = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar);
    
                % Create condition to find data
                depthCondition = currStatMonthData.depth >= zLower & ...
                                 currStatMonthData.depth <= zUpper;
            
                % CASE 1: Data found in window
                if sum(depthCondition) ~= 0
                    currStatMonthDhData = currStatMonthData(depthCondition,:);
                    depthMatched = unique(currStatMonthDhData.depth);
    
                    % If in the surface ocean, redefine zeu bounds based on matched depths
                    if iDh == 1 
                        if numel(depthMatched) == 1
                            % Only one depth — expand symmetrically ±10%
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(depthMatched - OC_ERR_FRAC * depthMatched);
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(depthMatched + OC_ERR_FRAC * depthMatched);

                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar);
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar);
                        elseif numel(depthMatched) > 1 
                            % Multiple depths — adapt based on relation to Zeu
                            zeu = qZeuMonthly(iMonth,iLoc);
                            if all(depthMatched < zeu)
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(min(depthMatched));
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(zeu);
                            elseif all(depthMatched > zeu)
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(zeu);
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(max(depthMatched));
                            else
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(min(depthMatched));
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(max(depthMatched));
                            end
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar);
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar);
                        end
        
                        % Update variables
                        zLower = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar);
                        zUpper = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar);
        
                        depthCondition = currStatMonthData.depth >= zLower & ...
                                         currStatMonthData.depth <= zUpper;
                        currStatMonthDhData = currStatMonthData(depthCondition,:);
                    end
    
                % CASE 2: No data found within the depth window
                elseif sum(depthCondition) == 0
                    
                    % If in the surface ocean, find closest available depth 
                    % data (closest depth that is < 200 m)
                    if iDh == 1 && any(currStatMonthData.depth < MAX_ZEU)
    
                        % Compute distance to Zeu
                        zeu = qZeuMonthly(iMonth,iLoc);
                        absDiff = abs(currStatMonthData.depth - zeu);
    
                        % First try a ±20m margin
                        currStatMonthDhData = currStatMonthData(absDiff < 20,:);
    
                        % If such data do not exist (i.e., the margin is
                        % larger than 20 m) extract data that match the 
                        % minimum absolute difference
                        if isempty(currStatMonthDhData)
                            minDiff = min(absDiff); % find the minimum absolute difference
                            currStatMonthDhData = currStatMonthData(absDiff == minDiff,:); % extract data that match the minimum absolute difference
                        end
                    
                        depthMatched = unique(currStatMonthDhData.depth);
    
                        % Redefine bounds based on matched depths
                        if numel(depthMatched) == 1 
                            % Only one depth — expand symmetrically ±10%
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(depthMatched - OC_ERR_FRAC * depthMatched);
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(depthMatched + OC_ERR_FRAC * depthMatched);

                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar);
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar);
                        elseif numel(depthMatched) > 1
                            % Multiple depths — adapt based on relation to Zeu
                            zeu = qZeuMonthly(iMonth,iLoc);
                            if all(depthMatched < zeu)
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(min(depthMatched));
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(zeu);
                            elseif all(depthMatched > zeu)
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(zeu);
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(max(depthMatched));
                            else
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = floor(min(depthMatched));
                                LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = ceil(max(depthMatched));
                            end
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar);
                            LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh,iVar) = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar);
                        end
    
                        % Update variables
                        zLower = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar);
                        zUpper = LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar);
        
                        depthCondition = currStatMonthData.depth >= zLower & ...
                                         currStatMonthData.depth <= zUpper;
                        currStatMonthDhData = currStatMonthData(depthCondition,:);
    
                    % No data at all for this depth/month/station (only update 
                    % LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT)
                    else
                        currStatMonthDhData = [];
                        LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) = NaN;
                        LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar) = NaN;
                    end
                end

                % Proceed if there are data in the current month
                if ~isempty(currStatMonthDhData)
    
                    % Add the depth horizon tag to currStatMonthDhData
                    currStatMonthDhData.depthHorizon(:) = tagDepthHorizons(iDh);
    
                    % Add the depth horizon tag to D array for later use
                    idxRows = find(D.tag == STATION_NAMES{iLoc} &...
                            D.depth >= LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,1,iDh,iVar) &...
                            D.depth <= LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(iMonth,iLoc,2,iDh,iVar));
                    D.depthHorizon(idxRows) = tagDepthHorizons(iDh);
    
                    % Store values, their depth and their type
                    picbsiFluxRawDhValues_cell{iDh,iMonth,iLoc,iVar} = num2str(currStatMonthDhData{:,1}'); 
                    picbsiFluxRawDhDepths_cell{iDh,iMonth,iLoc,iVar} = num2str(currStatMonthDhData{:,5}');
                    picbsiFluxRawDhDataType_cell{iDh,iMonth,iLoc,iVar} = strjoin(currStatMonthDhData{:,6}'); 
                    picbsiFluxRawDhTag_cell{iDh,iMonth,iLoc,iVar} = strjoin(currStatMonthDhData{:,7}');
    
                    % (4) Calculate the average & get the associated no. data points (N)
                    picbsiFluxMonthlyDhAvg(iDh,iMonth,iLoc,iVar) = mean(currStatMonthDhData{:,1}); % mmol m-2 d-1
                    picbsiFluxMonthlyDhN(iDh,iMonth,iLoc,iVar) = sum(currStatMonthDhData{:,1} >= 0,'omitnan');
                    
                    % (5) Calculate net error (propagate type A and type B errors separately)
                    vals = squeeze(currStatMonthDhData{:,1});
                    errSys = squeeze(currStatMonthDhData{:,2});
                    errRand = squeeze(currStatMonthDhData{:,3});
                    errSys(isnan(errSys)) = 0;
                    errRand(isnan(errRand)) = 0;
    
                    % Error propagation of random errors (type A)
                    [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),vals,errRand);
                    typeAuncertainty = f_UB - f_MID;
    
                    % Error propagation of the systematic errors (type B)
                    [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),vals,errSys);
                    typeBuncertainty = f_UB - f_MID;
    
                    % Net error: sum of type A and type B errors in the quadrature
                    picbsiFluxMonthlyDhErrTot(iDh,iMonth,iLoc,iVar) =...
                        sqrt(typeAuncertainty^2 + typeBuncertainty^2); % mmol m-2 d-1
    
                end % are there data in the current month?
            end % iDh    
        end % iMonth
    end % iLoc

    clear vals errRand errSys

end % iVar

% Tide up the variable 'depthHorizon'. Find the rows that have not been 
% assigned into a depth horizon and assign them a 'NaN' string
isEmptyDhTag = cellfun(@isempty,D.depthHorizon);
iEmptyRows = find(isEmptyDhTag);
nEmptyRows = length(iEmptyRows);
D.depthHorizon(iEmptyRows) = repmat({'NaN'},nEmptyRows,1);
D.depthHorizon = categorical(D.depthHorizon);
TRAPRAD_TABLE_PICBSI = D;
DP = D(D.depthHorizon ~= 'NaN',:); % data processed

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - BIN DATA MONTHLY BY UNIQUE DEPTH AND PROPAGATE ERROR
% -------------------------------------------------------------------------

% Unfold the data stored in three cell arrays and store them into numerical 
% arrays.

picbsiFluxRawProfileValues   = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS,NUM_VARS);
picbsiFluxRawProfileErrRand  = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS,NUM_VARS);
picbsiFluxRawProfileErrSys   = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS,NUM_VARS);
picbsiFluxRawProfileDepths   = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS,NUM_VARS);
picbsiFluxRawProfileDataType = cell(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS,NUM_VARS);

for iVar = 1:NUM_VARS
    for iLoc = 1:NUM_LOCS
        for iMonth = 1:12
    
            currMonthData     = zeros(MAX_NUM_VALUES_PER_MONTH,1);
            currMonthErrRand  = zeros(MAX_NUM_VALUES_PER_MONTH,1);
            currMonthErrSys   = zeros(MAX_NUM_VALUES_PER_MONTH,1);
            currMonthDepths   = zeros(MAX_NUM_VALUES_PER_MONTH,1);
            currMonthDataType = cell(MAX_NUM_VALUES_PER_MONTH,1);
    
            allmyvals = picbsiFluxRawProfileValues_cell{iMonth,iLoc,iVar};
            allmyranderrs = picbsiFluxRawProfileErrRand_cell{iMonth,iLoc,iVar};
            allmysyserrs = picbsiFluxRawProfileErrSys_cell{iMonth,iLoc,iVar};
            allmydepths = picbsiFluxRawProfileDepths_cell{iMonth,iLoc,iVar};
            allmydatatypes = picbsiFluxRawProfileDataType_cell{iMonth,iLoc,iVar};
            
            if (~isempty(allmyvals))
                nDataPoints = length(str2num(allmyvals));
                currMonthData(1:nDataPoints) = str2num(allmyvals);
                currMonthErrSys(1:nDataPoints) = str2num(allmysyserrs);
                currMonthErrRand(1:nDataPoints) = str2num(allmyranderrs);
                thetypesall = textscan(allmydatatypes,'%s');
                for iDataPoint = 1:nDataPoints
                    currMonthDataType{iDataPoint} = thetypesall{1}{iDataPoint};
                end
                currMonthDepths(1:nDataPoints) = str2num(allmydepths);
                
                % Correct for NaN
                theranderrs = currMonthErrRand(1:nDataPoints);
                theranderrs(isnan(theranderrs)) = 0;
                currMonthErrRand(1:nDataPoints) = theranderrs;
            else
                nDataPoints = 0;
            end
            
            [rowIdxs,~,yall] = find(currMonthData);
            xall = zeros(nDataPoints,1);
            gall = cell(nDataPoints,1);
            for iDataPoint = 1:nDataPoints
                xall(iDataPoint) = currMonthDepths(rowIdxs(iDataPoint));
                gall(iDataPoint) = currMonthDataType(rowIdxs(iDataPoint));
            end
            
            picbsiFluxRawProfileValues(1:nDataPoints,iMonth,iLoc,iVar) = yall; % mmol m-2 d-1
            picbsiFluxRawProfileErrSys(1:nDataPoints,iMonth,iLoc,iVar) = currMonthErrSys(1:nDataPoints); % mmol m-2 d-1
            picbsiFluxRawProfileErrRand(1:nDataPoints,iMonth,iLoc,iVar) = currMonthErrRand(1:nDataPoints); % mmol m-2 d-1
            picbsiFluxRawProfileDepths(1:nDataPoints,iMonth,iLoc,iVar) = xall;
            picbsiFluxRawProfileDataType(1:nDataPoints,iMonth,iLoc,iVar) = gall;
    
        end % iMonth 
    end % iLoc
end % iVar

% Average data per unique depth and propagate error.
 
picbsiFluxMonthlyProfileAvg    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS,NUM_VARS); % mean
picbsiFluxMonthlyProfileN      = zeros(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS,NUM_VARS); % number of values
picbsiFluxMonthlyProfileErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS,NUM_VARS); % total error
picbsiFluxMonthlyProfileDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS,NUM_VARS);  
nUniqueObsDepths               = NaN(12,NUM_LOCS,NUM_VARS); 

for iVar = 1:NUM_VARS
    for iLoc = 1:NUM_LOCS
        for iMonth = 1:12
    
            % Get values
            valso = squeeze(picbsiFluxRawProfileValues(:,iMonth,iLoc,iVar));
            errssyso = squeeze(picbsiFluxRawProfileErrSys(:,iMonth,iLoc,iVar));
            errsrando = squeeze(picbsiFluxRawProfileErrRand(:,iMonth,iLoc,iVar));
            zo = squeeze(picbsiFluxRawProfileDepths(:,iMonth,iLoc,iVar));
            zo(zo==0) = NaN;
            uniqueObsDepths = unique(zo);
            nUniqueObsDepths(iMonth,iLoc,iVar) = sum(~isnan(uniqueObsDepths));
            
            % Calculate average value and error per unique depth
            for iUniqueDepth = 1:nUniqueObsDepths(iMonth,iLoc,iVar)
                idxThisUniqueDepth = zo == uniqueObsDepths(iUniqueDepth);
                fluxesInThisUniqueDepth = valso(idxThisUniqueDepth);
                sysErrsInThisUniqueDepth = errssyso(idxThisUniqueDepth);
                randErrsInThisUniqueDepth = errsrando(idxThisUniqueDepth);
                
                % Depth value
                picbsiFluxMonthlyProfileDepths(iUniqueDepth,iMonth,iLoc,iVar) = uniqueObsDepths(iUniqueDepth);
                
                % Calculate the average and get the number of associated data points
                picbsiFluxMonthlyProfileAvg(iUniqueDepth,iMonth,iLoc,iVar) = mean(fluxesInThisUniqueDepth);
                picbsiFluxMonthlyProfileN(iUniqueDepth,iMonth,iLoc,iVar) = sum(fluxesInThisUniqueDepth >= 0,'omitnan');
                
                % Type A uncertainty: two options,
                % (a) if no errors associated to the samples, calculate as standard deviation/sqrt(N)
                if (sum(randErrsInThisUniqueDepth) == 0)
                    typeAuncertainty = std(fluxesInThisUniqueDepth)./sqrt(picbsiFluxMonthlyProfileN(iUniqueDepth,iMonth,iLoc,iVar)); 
                % (b) if error associated to the samples, propagate error
                else
                    [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),fluxesInThisUniqueDepth,randErrsInThisUniqueDepth);
                    typeAuncertainty = f_UB - f_MID;
                end
    
                % Type B uncertainty, or systematic error 
                % Error propagation of the systematic errors
                [~,~,~,f_MID,f_UB,~,~] = worstcase(@(x) mean(x),fluxesInThisUniqueDepth,sysErrsInThisUniqueDepth);
                typeBuncertainty = f_UB - f_MID;
    
                % Net error: sum of type A and type B errors in the quadrature
                picbsiFluxMonthlyProfileErrTot(iUniqueDepth,iMonth,iLoc,iVar) =... 
                    sqrt(typeAuncertainty^2 + typeBuncertainty^2);
      
            end % iUniqueDepth
        end % iMonth
    end % iLoc
end % iVar

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - BIN DATA ANNUALLY BY DEPTH HORIZON AND PROPAGATE ERROR
% -------------------------------------------------------------------------

picbsiFluxAnnualDhAvg    = NaN(NUM_TARGET_DEPTHS,NUM_LOCS,NUM_VARS); % weighted mean
picbsiFluxAnnualDhErrTot = NaN(NUM_TARGET_DEPTHS,NUM_LOCS,NUM_VARS);
picbsiFluxAnnualDhN      = NaN(NUM_TARGET_DEPTHS,NUM_LOCS,NUM_VARS);
picbsiFluxAnnualDhMin    = NaN(NUM_TARGET_DEPTHS,NUM_LOCS,NUM_VARS);
picbsiFluxAnnualDhMax    = NaN(NUM_TARGET_DEPTHS,NUM_LOCS,NUM_VARS);

for iVar = 1:NUM_VARS
    for iLoc = 1:NUM_LOCS
        for iDh = 1:NUM_TARGET_DEPTHS
            nSamples = squeeze(picbsiFluxMonthlyDhN(iDh,:,iLoc,iVar)); % nz x 12 x nLocs x nVars
            
            if any(nSamples)
                picbsiFluxAnnualDhN(iDh,iLoc,iVar) = sum(nSamples);
                
                vals = squeeze(picbsiFluxMonthlyDhAvg(iDh,:,iLoc,iVar));
                errTot = squeeze(picbsiFluxMonthlyDhErrTot(iDh,:,iLoc,iVar));
    
                % Set to 0 positions with no samples
                vals(nSamples == 0) = 0;
                errTot(nSamples == 0) = 0;
    
                % Calculate weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
                paramsWeightedAverage = vals.*nSamples;  
                calculateWeightedAvg = @(x) sum(x) ./ sum(nSamples(:), 'omitnan');
                picbsiFluxAnnualDhAvg(iDh,iLoc,iVar) = calculateWeightedAvg(paramsWeightedAverage);
    
                % Error propagation using worstcase with the function handle and the parameters
                [~,~,~,f_MID,f_UB,~,~] = ...
                    worstcase(@(paramsWeightedAverage) calculateWeightedAvg(paramsWeightedAverage),...
                    paramsWeightedAverage', errTot');
                picbsiFluxAnnualDhErrTot(iDh,iLoc,iVar) = f_UB - f_MID;
                
                % Min and max values
                picbsiFluxAnnualDhMin(iDh,iLoc,iVar) = min(picbsiFluxMonthlyDhAvg(iDh,:,iLoc,iVar));
                picbsiFluxAnnualDhMax(iDh,iLoc,iVar) = max(picbsiFluxMonthlyDhAvg(iDh,:,iLoc,iVar));
            end
        end
    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - BIN DATA ANNUALLY BY UNIQUE DEPTH AND PROPAGATE ERROR
% -------------------------------------------------------------------------

picbsiFluxAnnualProfileAvg    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,NUM_VARS); % weighted mean
picbsiFluxAnnualProfileErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,NUM_VARS);
picbsiFluxAnnualProfileN      = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,NUM_VARS);
picbsiFluxAnnualProfileMin    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,NUM_VARS);
picbsiFluxAnnualProfileMax    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,NUM_VARS);
picbsiFluxAnnualProfileDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS,NUM_VARS);  

for iVar = 1:NUM_VARS
    for iLoc = 1:NUM_LOCS
    
        % Get values
        vals = squeeze(picbsiFluxMonthlyProfileAvg(:,:,iLoc,iVar));
        errTot = squeeze(picbsiFluxMonthlyProfileErrTot(:,:,iLoc,iVar));
        nSamples = squeeze(picbsiFluxMonthlyProfileN(:,:,iLoc,iVar));
    
        zo = squeeze(picbsiFluxMonthlyProfileDepths(:,:,iLoc,iVar));
        zo(zo==0) = NaN;
        uniqueObsDepths = unique(zo);
        nUniqueObsDepths = sum(~isnan(uniqueObsDepths));
    
        % Calculate average value and error per unique depth
        for iUniqueDepth = 1:nUniqueObsDepths
    
            idxThisUniqueDepth = zo == uniqueObsDepths(iUniqueDepth);
            fluxesInThisUniqueDepth = vals(idxThisUniqueDepth);
            errsInThisUniqueDepth = errTot(idxThisUniqueDepth);
            nSamplesInThisUniqueDepth = nSamples(idxThisUniqueDepth);
            picbsiFluxAnnualProfileN(iUniqueDepth,iLoc,iVar) = sum(nSamplesInThisUniqueDepth);
    
            % Depth value
            picbsiFluxAnnualProfileDepths(iUniqueDepth,iLoc,iVar) = uniqueObsDepths(iUniqueDepth);
    
            % Calculate weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
            paramsWeightedAverage = fluxesInThisUniqueDepth.*nSamplesInThisUniqueDepth;  
            calculateWeightedAvg = @(x) sum(x) ./ sum(nSamplesInThisUniqueDepth(:), 'omitnan');
            picbsiFluxAnnualProfileAvg(iUniqueDepth,iLoc,iVar) = calculateWeightedAvg(paramsWeightedAverage);
    
            % Error propagation using worstcase with the function handle and the parameters
            [~,~,~,f_MID,f_UB,~,~] = ...
                worstcase(@(paramsWeightedAverage) calculateWeightedAvg(paramsWeightedAverage),...
                paramsWeightedAverage, errsInThisUniqueDepth);
            picbsiFluxAnnualProfileErrTot(iUniqueDepth,iLoc,iVar) = f_UB - f_MID;
    
            % Min and max values
            picbsiFluxAnnualProfileMin(iUniqueDepth,iLoc,iVar) = min(fluxesInThisUniqueDepth);
            picbsiFluxAnnualProfileMax(iUniqueDepth,iLoc,iVar) = max(fluxesInThisUniqueDepth);
    
        end % iUniqueDepth
    end % iLoc
end % iVar

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - CALCULATE THE NUMBER OF DATA POINTS BASED ON VARIOUS CRITERIA
% -------------------------------------------------------------------------

% Entries for depth horizon analysis
picbsiFluxMonthlyDhNbyMethod = zeros(12,NUM_LOCS,NUM_TARGET_DEPTHS,2,NUM_VARS); % same as classicMonthlyDhN but split by collection method
for iVar = 1:NUM_VARS
    for iLoc = 1:NUM_LOCS
        for iMonth = 1:12
            for iDh = 1:NUM_TARGET_DEPTHS
                myvals = picbsiFluxRawDhValues_cell{iDh,iMonth,iLoc,iVar};
                mytypes = picbsiFluxRawDhDataType_cell{iDh,iMonth,iLoc,iVar};
                if ~isempty(myvals)
                    thetypes = textscan(mytypes,'%s');
                    unnestthetypes = [thetypes{:}];
                    ntrap = nnz(strcmp(unnestthetypes,'trap'));
                    nradio = nnz(strcmp(unnestthetypes,'radionuclide'));
                    picbsiFluxMonthlyDhNbyMethod(iMonth,iLoc,iDh,1,iVar) = ntrap;
                    picbsiFluxMonthlyDhNbyMethod(iMonth,iLoc,iDh,2,iVar) = nradio;
                end
            end
        end
    end
end

for iVar = 1:NUM_VARS

    nObsForDhAnalysis = sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,:,:,:,iVar)),'all','omitnan'); % same as sum(squeeze(dataMonthlyN(:,:,:)),'all','omitnan');
    fprintf('\nWe have %d data points FOR SUMMARIES BY DEPTH HORIZON, of which', nObsForDhAnalysis)
    
    fracBatsDataPointsDhAnalysis  = (sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,4,:,:,iVar)),'all','omitnan')/nObsForDhAnalysis)*100;
    fracOspDataPointsDhAnalysis   = (sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,2,:,:,iVar)),'all','omitnan')/nObsForDhAnalysis)*100;
    fracPapsoDataPointsDhAnalysis = (sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,3,:,:,iVar)),'all','omitnan')/nObsForDhAnalysis)*100;
    fracEqpacDataPointsDhAnalysis = (sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,1,:,:,iVar)),'all','omitnan')/nObsForDhAnalysis)*100;
    fracHotalohaDataPointsDhAnalysis = (sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,5,:,:,iVar)),'all','omitnan')/nObsForDhAnalysis)*100;
    fracHausgartenDataPointsDhAnalysis = (sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,6,:,:,iVar)),'all','omitnan')/nObsForDhAnalysis)*100;
    
    fprintf('\n%0.1f%% data points at BATS/OFP,', fracBatsDataPointsDhAnalysis)
    fprintf('\n%0.1f%% data points at HOT/ALOHA, and', fracHotalohaDataPointsDhAnalysis)
    fprintf('\n%0.1f%% data points at OSP,', fracOspDataPointsDhAnalysis)
    fprintf('\n%0.1f%% data points at PAP-SO,', fracPapsoDataPointsDhAnalysis)
    fprintf('\n%0.1f%% data points at HAUSGARTEN.', fracHausgartenDataPointsDhAnalysis)
    fprintf('\n%0.1f%% data points at EqPac,', fracEqpacDataPointsDhAnalysis)
    
    nRadionuclideDataPoints = sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,:,:,2,iVar)),'all','omitnan');
    fracRadionuclideDataPoints = (nRadionuclideDataPoints/nObsForDhAnalysis)*100;
    fprintf('\nWe have %0.1f%% radionuclide data points for summaries by depth horizon.', fracRadionuclideDataPoints)

end % iVar

% Compare with length of the raw data set
nEntriesRawDatasetByLoc = zeros(NUM_LOCS,NUM_VARS);

for iVar = 1:NUM_VARS
    for iLoc = 1:NUM_LOCS
        if (iVar == 1)
            currStatData = D.PIC_mmol_m2_d(D.tag == STATION_NAMES{iLoc},:);
        else
            currStatData = D.Si_mmol_m2_d(D.tag == STATION_NAMES{iLoc},:);
        end
        nEntriesRawDatasetByLoc(iLoc,iVar) = nnz(~isnan(currStatData));
    end
    nEntriesRawDatasetTotal = sum(nEntriesRawDatasetByLoc(:,iVar)); % = nnz(~isnan(dataRaw.POC_mmol_m2_d));
    fprintf('\nWe have %d data points IN TOTAL, of which', nEntriesRawDatasetTotal)
    
    fracBatsDataPointsRaw   = (nEntriesRawDatasetByLoc(4,iVar)/nEntriesRawDatasetTotal)*100;
    fracOspDataPointsRaw    = (nEntriesRawDatasetByLoc(2,iVar)/nEntriesRawDatasetTotal)*100;
    fracPapsoDataPointsRaw  = (nEntriesRawDatasetByLoc(3,iVar)/nEntriesRawDatasetTotal)*100;
    fracEqpacDataPointsRaw  = (nEntriesRawDatasetByLoc(1,iVar)/nEntriesRawDatasetTotal)*100;
    fracHotalohaDataPointsRaw  = (nEntriesRawDatasetByLoc(5,iVar)/nEntriesRawDatasetTotal)*100;
    fracHusgartenDataPointsRaw  = (nEntriesRawDatasetByLoc(6,iVar)/nEntriesRawDatasetTotal)*100;
    
    fprintf('\n%0.1f%% data points at BATS/OFP', fracBatsDataPointsRaw)
    fprintf('\n%0.1f%% data points at HOT/ALOHA', fracHotalohaDataPointsRaw)
    fprintf('\n%0.1f%% data points at OSP', fracOspDataPointsRaw)
    fprintf('\n%0.1f%% data points at PAP-SO', fracPapsoDataPointsRaw)
    fprintf('\n%0.1f%% data points at HAUSGARTEN', fracHusgartenDataPointsRaw)
    fprintf('\n%0.1f%% data points at EqPac', fracEqpacDataPointsRaw)
    
    if (iVar == 1)
        nEntriesRawDatasetRadionuclideTotal = length(D.PIC_mmol_m2_d(strcmp(D.method,'radionuclide'),:));
    else
        nEntriesRawDatasetRadionuclideTotal = length(D.Si_mmol_m2_d(strcmp(D.method,'radionuclide'),:));  
    end
    fracRadionuclidesDatasetRaw = (nEntriesRawDatasetRadionuclideTotal/nEntriesRawDatasetTotal)*100;
    fprintf('\nWe have %0.1f%% radionuclide data points in total.\n', fracRadionuclidesDatasetRaw)

end

% Print number of data points for reference
tablendp = zeros(NUM_TARGET_DEPTHS,NUM_LOCS,2,NUM_VARS);
for iVar = 1:NUM_VARS
    for iLoc = 1:NUM_LOCS
        for iDh = 1:NUM_TARGET_DEPTHS
            for iMethod = 1:2
                ndp = sum(squeeze(picbsiFluxMonthlyDhNbyMethod(:,iLoc,iDh,iMethod,iVar)),'omitnan');
                tablendp(iDh,iLoc,iMethod,iVar) = ndp;
            end
        end
    end
end
checkNumEntries = tablendp(1,6,2,1); % see in 6th loc, zeu for method 2 (radionuclide), 1=PIC

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 8 - SAVE THE DATA
% -------------------------------------------------------------------------

save(fullfile('.','data','processed',filenameOutputFluxCompilation),...
    'picbsiFluxRawProfileValues_cell','picbsiFluxRawProfileErrRand_cell',...
    'picbsiFluxRawProfileErrSys_cell','picbsiFluxRawProfileDepths_cell',...
    'picbsiFluxRawProfileDataType_cell','picbsiFluxRawProfileValues',...
    'picbsiFluxRawProfileErrRand','picbsiFluxRawProfileErrSys',...
    'picbsiFluxRawProfileDepths','picbsiFluxRawProfileDataType',...
    'picbsiFluxRawDhValues_cell','picbsiFluxRawDhDepths_cell',...
    'picbsiFluxRawDhTag_cell','picbsiFluxRawDhDataType_cell',...
    'picbsiFluxMonthlyProfileAvg','picbsiFluxMonthlyProfileErrTot',...
    'picbsiFluxMonthlyProfileN','picbsiFluxMonthlyProfileDepths',...
    'picbsiFluxMonthlyDhAvg','picbsiFluxMonthlyDhErrTot','picbsiFluxMonthlyDhN',...
    'picbsiFluxAnnualDhAvg','picbsiFluxAnnualDhErrTot','picbsiFluxAnnualDhN',...
    'picbsiFluxAnnualProfileAvg','picbsiFluxAnnualProfileErrTot',...
    'picbsiFluxAnnualProfileN','picbsiFluxAnnualProfileDepths',...
    'TRAPRAD_TABLE_PICBSI','LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT',...
    'LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT_FILLED','-v7.3')

fprintf('\nThe PIC and bSi flux compilation data have been saved correctly.\n')