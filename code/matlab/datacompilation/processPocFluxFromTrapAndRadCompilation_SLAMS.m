
% ======================================================================= %
%                                                                         %
% This script reads in the POC flux compilation created for the GRL 2025  %
% publication (Dataset S0) but modified to suit calibration purposes of   %
% SLAMS-2.0, where 4 depth horizons of data are needed (instead of the    %
% default 3). It then calculates monthly and annual averages and          %
% propagates error accordingly. Output units are "mmol C m-2 d-1".        %
%                                                                         %
% The script has 9 sections:                                              %
%   Section 1 - Presets.                                                  %
%   Section 2 - Load the dataset and manipulate the data array.           %
%   Section 3 - Get euphotic layer depth.                                 %
%   Section 4 - Bin data monthly by depth horizon and propagate error.    %
%   Section 5 - Bin data monthly by unique depth and propagate error.     %
%   Section 6 - Bin data annually by depth horizon and propagate error.   %
%   Section 7 - Bin data annually by unique depth and propagate error.    %
%   Section 8 - Calculate the number of data points based on various      % 
%               criteria.                                                 %
%   Section 9 - Save the data.                                            %              
%                                                                         %
%   This script uses these external functions:                            % 
%       fillNansStrategicallyInSurfaceOceanData.m - custom function       %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 8 May 2025                                    %                                  
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/raw/'));
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
filenameInputFluxCompilation        = 'dataset_s0_trap_and_radionuclide_compilation.xlsx';
filenameInputZeu                    = 'zeu_calculated_chlaoccci_mldifremer_pointonepercentpar0.mat';
filenameOutputFluxCompilation       = 'pocflux_compilation_slams.mat';
filenameOutputTimeseriesInformation = 'timeseries_station_information_slams.mat';

% Log progress
logID = fopen(fullfile('.','data','interim','log_pocflux_compilation.log'),'w'); 

% Define parameters
POCFLUX_RAND_ERR_FRAC = 0.30; % 30% (Buesseler et al. 2000, Buesseler et al. 2007, Stanley et al. 2004)
POCFLUX_SYS_ERR_FRAC = 0.10; % 10%, based on a literature review
OC_ERR_FRAC = 0.10; % ocean colour error fraction (10%), based on McKinna et al. 2019
MAX_ZEU = 200; % m
MOLAR_MASS_CARBON = 12.011; % g mol-1
MAX_NUM_VALUES_PER_MONTH = 1000;
MAX_NUM_DEPTHS_PER_PROFILE = 100;

% Load the global-ocean euphotic layer depth product
load(fullfile('.','data','interim',filenameInputZeu),'zeu','zeu_lat','zeu_lon')

% Station information
STATION_NAMES = {'EqPac','OSP','PAP-SO','BATS/OFP','HOT/ALOHA','HAUSGARTEN'}; 
STATION_TAGS = {'eqpac','osp','papso','batsofp','hotaloha','hausgarten'}; 
NUM_LOCS = length(STATION_NAMES);

% The raw POC flux data array should contain data for the following ocean
% sites, which have the following latitude and longitudes:

          % EqPac  OSP PAP-SO   BATS  HOT HAUSGARTEN
LOC_LATS = [   0,   50,    49,  31.6, 22.5,  79]; % HOT/ALOHA lat slightly modified to extract data from WOA 
LOC_LONS = [-140, -145, -16.5, -64.2, -158, 4.3];

% The following depth horizons defined on the data will be used for data 
% extraction
tagDepthHorizons = {'zeu','mesoupp','mesolow','zmeso'};
NUM_TARGET_DEPTHS = length(tagDepthHorizons);
LOC_DEPTH_HORIZONS = zeros(12,NUM_LOCS,2,NUM_TARGET_DEPTHS); % 2 boundaries for each of the 4 depth horizons

% EqPac                              % OSP                                % PAP-SO                          
% Targets data at zeu                Targets data at zeu                  Targets data at zeu
LOC_DEPTH_HORIZONS(:,1,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,2,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,3,1,1) = NaN;   
LOC_DEPTH_HORIZONS(:,1,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,2,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,3,2,1) = NaN;   

% Targets the traps at 150-200 m     Targets the traps at 250-300 m       Targets the traps at 171-200 m     
LOC_DEPTH_HORIZONS(:,1,1,2) = 150;   LOC_DEPTH_HORIZONS(:,2,1,2) = 250;   LOC_DEPTH_HORIZONS(:,3,1,2) = 171;  
LOC_DEPTH_HORIZONS(:,1,2,2) = 200;   LOC_DEPTH_HORIZONS(:,2,2,2) = 300;   LOC_DEPTH_HORIZONS(:,3,2,2) = 200; 

% Targets the trap at 320 m          Targets the trap at 600 m            Targets the traps at 446-500 m     
LOC_DEPTH_HORIZONS(:,1,1,3) = 310;   LOC_DEPTH_HORIZONS(:,2,1,3) = 590;   LOC_DEPTH_HORIZONS(:,3,1,3) = 446;  
LOC_DEPTH_HORIZONS(:,1,2,3) = 330;   LOC_DEPTH_HORIZONS(:,2,2,3) = 610;   LOC_DEPTH_HORIZONS(:,3,2,3) = 500; 

% Targets the trap at 880 m          Targets the traps at 1000 & 1009 m   Targets the trap at 1000 m     
LOC_DEPTH_HORIZONS(:,1,1,4) = 870;   LOC_DEPTH_HORIZONS(:,2,1,4) = 990;   LOC_DEPTH_HORIZONS(:,3,1,4) = 990;  
LOC_DEPTH_HORIZONS(:,1,2,4) = 890;   LOC_DEPTH_HORIZONS(:,2,2,4) = 1010;  LOC_DEPTH_HORIZONS(:,3,2,4) = 1010; 


% BATS/OFP                           % HOT/ALOHA                          % HAUSGARTEN
% Targets data at zeu                Targets data at zeu                  Targets data at zeu 
LOC_DEPTH_HORIZONS(:,4,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,5,1,1) = NaN;   LOC_DEPTH_HORIZONS(:,6,1,1) = NaN;
LOC_DEPTH_HORIZONS(:,4,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,5,2,1) = NaN;   LOC_DEPTH_HORIZONS(:,6,2,1) = NaN;

% Targets the trap at 300 m          Targets the trap at 300 m            Targets the traps at 280-314 m
LOC_DEPTH_HORIZONS(:,4,1,2) = 290;   LOC_DEPTH_HORIZONS(:,5,1,2) = 290;   LOC_DEPTH_HORIZONS(:,6,1,2) = 270;
LOC_DEPTH_HORIZONS(:,4,2,2) = 310;   LOC_DEPTH_HORIZONS(:,5,2,2) = 310;   LOC_DEPTH_HORIZONS(:,6,2,2) = 320;

% Targets the trap at 500 m          Targets the trap at 800 m            Targets the trap at 800 m 
LOC_DEPTH_HORIZONS(:,4,1,3) = 490;   LOC_DEPTH_HORIZONS(:,5,1,3) = 790;   LOC_DEPTH_HORIZONS(:,6,1,3) = 790;
LOC_DEPTH_HORIZONS(:,4,2,3) = 510;   LOC_DEPTH_HORIZONS(:,5,2,3) = 810;   LOC_DEPTH_HORIZONS(:,6,2,3) = 810;

% Targets the trap at 1500 m         Targets the trap at 1500 m           Targets the traps at 1225 & 1250 m
LOC_DEPTH_HORIZONS(:,4,1,4) = 1490;  LOC_DEPTH_HORIZONS(:,5,1,4) = 1490;  LOC_DEPTH_HORIZONS(:,6,1,4) = 1220;
LOC_DEPTH_HORIZONS(:,4,2,4) = 1510;  LOC_DEPTH_HORIZONS(:,5,2,4) = 1510;  LOC_DEPTH_HORIZONS(:,6,2,4) = 1260;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD THE POC FLUX DATASET AND MANIPULATE THE DATA ARRAY
% -------------------------------------------------------------------------

% Load the excel spreadhseet with POC flux data
opts = detectImportOptions(fullfile('..','Rufas_etal_2024_GRL','data','raw',filenameInputFluxCompilation));
opts = setvartype(opts,{'POC_mmol_m2_d','POC_mg_m2_d','randerr_POC_mmol_m2_d','randerr_POC_mg_m2_d'},'double');
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
D.syserr_POC_mmol_m2_d = POCFLUX_SYS_ERR_FRAC .* D.POC_mmol_m2_d;
D.syserr_POC_mg_m2_d = MOLAR_MASS_CARBON .* D.syserr_POC_mmol_m2_d;

% Update the random error
D.randerr_POC_mmol_m2_d(isnan(D.randerr_POC_mmol_m2_d)) = POCFLUX_RAND_ERR_FRAC .* D.POC_mmol_m2_d(isnan(D.randerr_POC_mmol_m2_d));
D.randerr_POC_mg_m2_d = MOLAR_MASS_CARBON .* D.randerr_POC_mmol_m2_d;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - GET EUPHOTIC LAYER DEPTH
% -------------------------------------------------------------------------

% Query points for interpolation
qLats = LOC_LATS;
qLons = LOC_LONS;

% Original data grid
[X,Y,T] = ndgrid(zeu_lat,zeu_lon,(1:12)');

% Interpolant 
F = griddedInterpolant(X, Y, T, zeu, 'linear'); 

% Extract data for the study locations (defined by qLats and qLons)
qZeuMonthly = NaN(12,NUM_LOCS);
for iLoc = 1:NUM_LOCS
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
    
    % Add information to LOC_DEPTH_HORIZONS array
    LOC_DEPTH_HORIZONS(:,iLoc,1,1) = floor(qZeuMonthly(:,iLoc)... 
        - (OC_ERR_FRAC.*qZeuMonthly(:,iLoc))); % -10%
    LOC_DEPTH_HORIZONS(:,iLoc,2,1) = ceil(qZeuMonthly(:,iLoc)... 
        + (OC_ERR_FRAC.*qZeuMonthly(:,iLoc))); % +10%
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - BIN DATA MONTHLY BY DEPTH HORIZON AND PROPAGATE ERROR (THESE
% DATA WILL BE USED FOR SUMMARIES AND PLOTTING)
% -------------------------------------------------------------------------

% Copy to be later on modified
LOC_DEPTH_HORIZONS_POC_OBS_ADAPT = LOC_DEPTH_HORIZONS;
LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED = LOC_DEPTH_HORIZONS;

% The variables that we want to extract
tagPocValues = 'POC_mmol_m2_d';
tagSysError = 'syserr_POC_mmol_m2_d';
tagRandError = 'randerr_POC_mmol_m2_d';

% Define output arrays
pocFluxRawProfileValues_cell   = cell(12,NUM_LOCS); 
pocFluxRawProfileErrRand_cell  = cell(12,NUM_LOCS); 
pocFluxRawProfileErrSys_cell   = cell(12,NUM_LOCS);
pocFluxRawProfileDepths_cell   = cell(12,NUM_LOCS); 
pocFluxRawProfileDataType_cell = cell(12,NUM_LOCS); 

pocFluxRawDhValues_cell   = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS); 
pocFluxRawDhDepths_cell   = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS); 
pocFluxRawDhDataType_cell = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS); % sediment trap vs radionuclides
pocFluxRawDhTag_cell      = cell(NUM_TARGET_DEPTHS,12,NUM_LOCS); % zeu, mesoupp, mesolow or zmeso

pocFluxMonthlyDhAvg    = NaN(NUM_TARGET_DEPTHS,12,NUM_LOCS);   % mean
pocFluxMonthlyDhN      = zeros(NUM_TARGET_DEPTHS,12,NUM_LOCS); % number of values
pocFluxMonthlyDhErrTot = NaN(NUM_TARGET_DEPTHS,12,NUM_LOCS);   % total error

for iLoc = 1:NUM_LOCS
    currStatData = D(D.tag == STATION_NAMES{iLoc},:);

    % Locate the variables that are relevant
    colIdsAll = zeros(6,1);
    colIdsAll(1) = find(strcmpi(currStatData.Properties.VariableNames,tagPocValues));
    colIdsAll(2) = find(strcmpi(currStatData.Properties.VariableNames,tagSysError));
    colIdsAll(3) = find(strcmpi(currStatData.Properties.VariableNames,tagRandError));
    colIdsAll(4) = find(strcmpi(currStatData.Properties.VariableNames,'month'));
    colIdsAll(5) = find(strcmpi(currStatData.Properties.VariableNames,'depth'));
    colIdsAll(6) = find(strcmpi(currStatData.Properties.VariableNames,'method'));

    % Extract these variables
    currStatData = currStatData(:,colIdsAll');
            
    % Identify rows with any NaN values
    rowsWithNaN = any(ismissing(currStatData), 2);

    % Remove rows with NaN values
    currStatDataValid = currStatData(~rowsWithNaN, :);
    
    % (1) Pull all the data at that station and save it to output arrays
    for iMonth = 1:12
        currMonthData = currStatDataValid(currStatDataValid.month == num2str(iMonth),:);
        if (~isempty(currMonthData))
            pocFluxRawProfileValues_cell{iMonth,iLoc} = num2str(currMonthData{:,1}'); 
            pocFluxRawProfileErrSys_cell{iMonth,iLoc} = num2str(currMonthData{:,2}');
            pocFluxRawProfileErrRand_cell{iMonth,iLoc} = num2str(currMonthData{:,3}'); 
            pocFluxRawProfileDepths_cell{iMonth,iLoc} = num2str(currMonthData{:,5}');
            pocFluxRawProfileDataType_cell{iMonth,iLoc} = strjoin(currMonthData{:,6}'); 
        end  
    end

    % (2) Pull data by month
    for iMonth = 1:12
        monthCondition = currStatDataValid.month == num2str(iMonth);
        currStatMonthData = currStatDataValid(monthCondition,:);
            
        % (3) Pull data by the depth horizon of interest
        for iDh = 1:NUM_TARGET_DEPTHS

            % Extract original depth bounds for the layer
            zLower = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh);
            zUpper = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh);

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
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(depthMatched - OC_ERR_FRAC * depthMatched);
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(depthMatched + OC_ERR_FRAC * depthMatched);
                        
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh);
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh);
                    elseif numel(depthMatched) > 1
                        % Multiple depths — adapt based on relation to Zeu
                        zeu = qZeuMonthly(iMonth,iLoc);
                        if all(depthMatched < zeu)
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(min(depthMatched));
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(zeu);
                        elseif all(depthMatched > zeu)
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(zeu);
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(max(depthMatched));
                        else
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(min(depthMatched));
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(max(depthMatched));
                        end
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh);
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh);
                    end

                    % Update variables
                    zLower = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh);
                    zUpper = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh);
    
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
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(depthMatched - OC_ERR_FRAC * depthMatched);
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(depthMatched + OC_ERR_FRAC * depthMatched);

                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh);
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh);
                    elseif numel(depthMatched) > 1
                        % Multiple depths — adapt based on relation to Zeu
                        zeu = qZeuMonthly(iMonth,iLoc);
                        if all(depthMatched < zeu)
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(min(depthMatched));
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(zeu);
                        elseif all(depthMatched > zeu)
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(zeu);
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(max(depthMatched));
                        else
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = floor(min(depthMatched));
                            LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = ceil(max(depthMatched));
                        end

                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,1,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh);
                        LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(iMonth,iLoc,2,iDh) = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh);
                    end

                    % Update variables
                    zLower = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh);
                    zUpper = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh);
    
                    depthCondition = currStatMonthData.depth >= zLower & ...
                                     currStatMonthData.depth <= zUpper;
                    currStatMonthDhData = currStatMonthData(depthCondition,:);

                % No data at all for this depth/month/station (only update 
                % LOC_DEPTH_HORIZONS_POC_OBS_ADAPT)
                else
                    currStatMonthDhData = [];
                    LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) = NaN;
                    LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh) = NaN;
                end
            end
            
            % Proceed if there are data in the current month
            if ~isempty(currStatMonthDhData)

                % Add the depth horizon tag to currStatMonthDhData
                currStatMonthDhData.depthHorizon(:) = tagDepthHorizons(iDh);

                % Add the depth horizon tag to D array for later use
                idxRows = find(D.tag == STATION_NAMES{iLoc} &...
                        D.depth >= LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,1,iDh) &...
                        D.depth <= LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(iMonth,iLoc,2,iDh));
                D.depthHorizon(idxRows) = tagDepthHorizons(iDh);

                % Store values, their depth and their type
                pocFluxRawDhValues_cell{iDh,iMonth,iLoc} = num2str(currStatMonthDhData{:,1}'); 
                pocFluxRawDhDepths_cell{iDh,iMonth,iLoc} = num2str(currStatMonthDhData{:,5}');
                pocFluxRawDhDataType_cell{iDh,iMonth,iLoc} = strjoin(currStatMonthDhData{:,6}'); 
                pocFluxRawDhTag_cell{iDh,iMonth,iLoc} = strjoin(currStatMonthDhData{:,7}');

                % (4) Calculate the average & get the associated no. data points (N)
                pocFluxMonthlyDhAvg(iDh,iMonth,iLoc) = mean(currStatMonthDhData{:,1}); % mmol C m-2 d-1
                pocFluxMonthlyDhN(iDh,iMonth,iLoc) = sum(currStatMonthDhData{:,1} >= 0,'omitnan');
                
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
                pocFluxMonthlyDhErrTot(iDh,iMonth,iLoc) =...
                    sqrt(typeAuncertainty^2 + typeBuncertainty^2); % mmol C m-2 d-1

            end % are there data in the current month?
        end % iDh    
    end % iMonth
end % iLoc

% % Checks
% a_vals = squeeze(classicMonthlyDhAvg(1,:,:));
% a_err = squeeze(classicMonthlyDhErrTot(1,:,:));
% a_stdperc = (a_err./a_vals).*100;
% a_n = squeeze(classicMonthlyDhN(NUM_LOCS,:,:));

% Tide up the variable 'depthHorizon'. Find the rows that have not been 
% assigned into a depth horizon and assign them a 'NaN' string
isEmptyDhTag = cellfun(@isempty,D.depthHorizon);
iEmptyRows = find(isEmptyDhTag);
nEmptyRows = length(iEmptyRows);
D.depthHorizon(iEmptyRows) = repmat({'NaN'},nEmptyRows,1);
D.depthHorizon = categorical(D.depthHorizon);
TRAPRAD_TABLE_POC = D;
DP = D(D.depthHorizon ~= 'NaN',:); % data processed

clear vals errRand errSys

% Save for use in other scripts
save(fullfile('.','data','interim',filenameOutputTimeseriesInformation),...
    'LOC_LATS','LOC_LONS','LOC_DEPTH_HORIZONS',...
    'STATION_NAMES','STATION_TAGS','MAX_NUM_VALUES_PER_MONTH',...
    'MAX_NUM_DEPTHS_PER_PROFILE','qZeuMonthly','MAX_ZEU')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - BIN DATA MONTHLY BY UNIQUE DEPTH AND PROPAGATE ERROR
% -------------------------------------------------------------------------

% Unfold the data stored in three cell arrays and store them into numerical 
% arrays.

pocFluxRawProfileValues   = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS);
pocFluxRawProfileErrRand  = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS);
pocFluxRawProfileErrSys   = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS);
pocFluxRawProfileDepths   = zeros(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS);
pocFluxRawProfileDataType = cell(MAX_NUM_VALUES_PER_MONTH,12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
    for iMonth = 1:12

        currMonthData     = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthErrRand  = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthErrSys   = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthDepths   = zeros(MAX_NUM_VALUES_PER_MONTH,1);
        currMonthDataType = cell(MAX_NUM_VALUES_PER_MONTH,1);

        allmyvals = pocFluxRawProfileValues_cell{iMonth,iLoc};
        allmyranderrs = pocFluxRawProfileErrRand_cell{iMonth,iLoc};
        allmysyserrs = pocFluxRawProfileErrSys_cell{iMonth,iLoc};
        allmydepths = pocFluxRawProfileDepths_cell{iMonth,iLoc};
        allmydatatypes = pocFluxRawProfileDataType_cell{iMonth,iLoc};
        
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
        
        pocFluxRawProfileValues(1:nDataPoints,iMonth,iLoc) = yall; % mmol C m-2 d-1
        pocFluxRawProfileErrSys(1:nDataPoints,iMonth,iLoc) = currMonthErrSys(1:nDataPoints); % mmol C m-2 d-1
        pocFluxRawProfileErrRand(1:nDataPoints,iMonth,iLoc) = currMonthErrRand(1:nDataPoints); % mmol C m-2 d-1
        pocFluxRawProfileDepths(1:nDataPoints,iMonth,iLoc) = xall;
        pocFluxRawProfileDataType(1:nDataPoints,iMonth,iLoc) = gall;

    end % iMonth 
end % iLoc

% Average data per unique depth and propagate error.
 
pocFluxMonthlyProfileAvg    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS); % mean
pocFluxMonthlyProfileN      = zeros(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS); % number of values
pocFluxMonthlyProfileErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS); % total error
pocFluxMonthlyProfileDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,12,NUM_LOCS);  
nUniqueObsDepths            = NaN(12,NUM_LOCS); 

for iLoc = 1:NUM_LOCS
    for iMonth = 1:12

        % Get values
        valso = squeeze(pocFluxRawProfileValues(:,iMonth,iLoc));
        errssyso = squeeze(pocFluxRawProfileErrSys(:,iMonth,iLoc));
        errsrando = squeeze(pocFluxRawProfileErrRand(:,iMonth,iLoc));
        zo = squeeze(pocFluxRawProfileDepths(:,iMonth,iLoc));
        zo(zo==0) = NaN;
        uniqueObsDepths = unique(zo);
        nUniqueObsDepths(iMonth,iLoc) = sum(~isnan(uniqueObsDepths));
        
        % Calculate average value and error per unique depth
        for iUniqueDepth = 1:nUniqueObsDepths(iMonth,iLoc)
            idxThisUniqueDepth = zo == uniqueObsDepths(iUniqueDepth);
            fluxesInThisUniqueDepth = valso(idxThisUniqueDepth);
            sysErrsInThisUniqueDepth = errssyso(idxThisUniqueDepth);
            randErrsInThisUniqueDepth = errsrando(idxThisUniqueDepth);
            
            % Depth value
            pocFluxMonthlyProfileDepths(iUniqueDepth,iMonth,iLoc) = uniqueObsDepths(iUniqueDepth);
            
            % Calculate the average and get the number of associated data points
            pocFluxMonthlyProfileAvg(iUniqueDepth,iMonth,iLoc) = mean(fluxesInThisUniqueDepth);
            pocFluxMonthlyProfileN(iUniqueDepth,iMonth,iLoc) = sum(fluxesInThisUniqueDepth >= 0,'omitnan');
            
            % Type A uncertainty: two options,
            % (a) if no errors associated to the samples, calculate as standard deviation/sqrt(N)
            if (sum(randErrsInThisUniqueDepth) == 0)
                typeAuncertainty = std(fluxesInThisUniqueDepth)./sqrt(pocFluxMonthlyProfileN(iUniqueDepth,iMonth,iLoc)); 
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
            pocFluxMonthlyProfileErrTot(iUniqueDepth,iMonth,iLoc) =... 
                sqrt(typeAuncertainty^2 + typeBuncertainty^2);
  
        end % iUniqueDepth
    end % iMonth
end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - BIN DATA ANNUALLY BY DEPTH HORIZON AND PROPAGATE ERROR
% -------------------------------------------------------------------------

pocFluxAnnualDhAvg    = NaN(NUM_TARGET_DEPTHS,NUM_LOCS); % weighted mean
pocFluxAnnualDhErrTot = NaN(NUM_TARGET_DEPTHS,NUM_LOCS);
pocFluxAnnualDhN      = NaN(NUM_TARGET_DEPTHS,NUM_LOCS);
pocFluxAnnualDhMin    = NaN(NUM_TARGET_DEPTHS,NUM_LOCS);
pocFluxAnnualDhMax    = NaN(NUM_TARGET_DEPTHS,NUM_LOCS);

for iLoc = 1:NUM_LOCS
    for iDh = 1:NUM_TARGET_DEPTHS
        nSamples = squeeze(pocFluxMonthlyDhN(iDh,:,iLoc)); % nz x 12 x nLocs
        
        if any(nSamples)
            pocFluxAnnualDhN(iDh,iLoc) = sum(nSamples);
            
            vals = squeeze(pocFluxMonthlyDhAvg(iDh,:,iLoc));
            errTot = squeeze(pocFluxMonthlyDhErrTot(iDh,:,iLoc));

            % Set to 0 positions with no samples
            vals(nSamples == 0) = 0;
            errTot(nSamples == 0) = 0;

            % Calculate weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
            paramsWeightedAverage = vals.*nSamples;  
            calculateWeightedAvg = @(x) sum(x) ./ sum(nSamples(:), 'omitnan');
            pocFluxAnnualDhAvg(iDh,iLoc) = calculateWeightedAvg(paramsWeightedAverage);

            % Error propagation using worstcase with the function handle and the parameters
            [~,~,~,f_MID,f_UB,~,~] = ...
                worstcase(@(paramsWeightedAverage) calculateWeightedAvg(paramsWeightedAverage),...
                paramsWeightedAverage', errTot');
            pocFluxAnnualDhErrTot(iDh,iLoc) = f_UB - f_MID;
            
            % Min and max values
            pocFluxAnnualDhMin(iDh,iLoc) = min(pocFluxMonthlyDhAvg(iDh,:,iLoc));
            pocFluxAnnualDhMax(iDh,iLoc) = max(pocFluxMonthlyDhAvg(iDh,:,iLoc));
        end
    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - BIN DATA ANNUALLY BY UNIQUE DEPTH AND PROPAGATE ERROR
% -------------------------------------------------------------------------

pocFluxAnnualProfileAvg    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS); % weighted mean
pocFluxAnnualProfileErrTot = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);
pocFluxAnnualProfileN      = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);
pocFluxAnnualProfileMin    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);
pocFluxAnnualProfileMax    = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);
pocFluxAnnualProfileDepths = NaN(MAX_NUM_DEPTHS_PER_PROFILE,NUM_LOCS);  

for iLoc = 1:NUM_LOCS

    % Get values
    vals = squeeze(pocFluxMonthlyProfileAvg(:,:,iLoc));
    errTot = squeeze(pocFluxMonthlyProfileErrTot(:,:,iLoc));
    nSamples = squeeze(pocFluxMonthlyProfileN(:,:,iLoc));

    zo = squeeze(pocFluxMonthlyProfileDepths(:,:,iLoc));
    zo(zo==0) = NaN;
    uniqueObsDepths = unique(zo);
    nUniqueObsDepths = sum(~isnan(uniqueObsDepths));

    % Calculate average value and error per unique depth
    for iUniqueDepth = 1:nUniqueObsDepths

        idxThisUniqueDepth = zo == uniqueObsDepths(iUniqueDepth);
        fluxesInThisUniqueDepth = vals(idxThisUniqueDepth);
        errsInThisUniqueDepth = errTot(idxThisUniqueDepth);
        nSamplesInThisUniqueDepth = nSamples(idxThisUniqueDepth);
        pocFluxAnnualProfileN(iUniqueDepth,iLoc) = sum(nSamplesInThisUniqueDepth);

        % Depth value
        pocFluxAnnualProfileDepths(iUniqueDepth,iLoc) = uniqueObsDepths(iUniqueDepth);

        % Calculate weighted mean (mw = ((mA*nA)+(mB*nB)+(mC*nC))/(nA+nB+nC))
        paramsWeightedAverage = fluxesInThisUniqueDepth.*nSamplesInThisUniqueDepth;  
        calculateWeightedAvg = @(x) sum(x) ./ sum(nSamplesInThisUniqueDepth(:), 'omitnan');
        pocFluxAnnualProfileAvg(iUniqueDepth,iLoc) = calculateWeightedAvg(paramsWeightedAverage);

        % Error propagation using worstcase with the function handle and the parameters
        [~,~,~,f_MID,f_UB,~,~] = ...
            worstcase(@(paramsWeightedAverage) calculateWeightedAvg(paramsWeightedAverage),...
            paramsWeightedAverage, errsInThisUniqueDepth);
        pocFluxAnnualProfileErrTot(iUniqueDepth,iLoc) = f_UB - f_MID;

        % Min and max values
        pocFluxAnnualProfileMin(iUniqueDepth,iLoc) = min(fluxesInThisUniqueDepth);
        pocFluxAnnualProfileMax(iUniqueDepth,iLoc) = max(fluxesInThisUniqueDepth);

    end % iUniqueDepth
end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 8 - CALCULATE THE NUMBER OF DATA POINTS BASED ON VARIOUS CRITERIA
% -------------------------------------------------------------------------

% Entries for depth horizon analysis
pocFluxMonthlyDhNbyMethod = zeros(12,NUM_LOCS,NUM_TARGET_DEPTHS,2); % same as classicMonthlyDhN but split by collection method
for iLoc = 1:NUM_LOCS
    for iMonth = 1:12
        for iDh = 1:NUM_TARGET_DEPTHS
            myvals = pocFluxRawDhValues_cell{iDh,iMonth,iLoc};
            mytypes = pocFluxRawDhDataType_cell{iDh,iMonth,iLoc};
            if ~isempty(myvals)
                thetypes = textscan(mytypes,'%s');
                unnestthetypes = [thetypes{:}];
                ntrap = nnz(strcmp(unnestthetypes,'trap'));
                nradio = nnz(strcmp(unnestthetypes,'radionuclide'));
                pocFluxMonthlyDhNbyMethod(iMonth,iLoc,iDh,1) = ntrap;
                pocFluxMonthlyDhNbyMethod(iMonth,iLoc,iDh,2) = nradio;
            end
        end
    end
end

nObsForDhAnalysis = sum(squeeze(pocFluxMonthlyDhNbyMethod(:,:,:,:)),'all','omitnan'); % same as sum(squeeze(dataMonthlyN(:,:,:)),'all','omitnan');
fprintf('\nWe have %d data points FOR SUMMARIES BY DEPTH HORIZON, of which', nObsForDhAnalysis)

fracBatsDataPointsDhAnalysis  = (sum(squeeze(pocFluxMonthlyDhNbyMethod(:,4,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracOspDataPointsDhAnalysis   = (sum(squeeze(pocFluxMonthlyDhNbyMethod(:,2,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracPapsoDataPointsDhAnalysis = (sum(squeeze(pocFluxMonthlyDhNbyMethod(:,3,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracEqpacDataPointsDhAnalysis = (sum(squeeze(pocFluxMonthlyDhNbyMethod(:,1,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracHotalohaDataPointsDhAnalysis = (sum(squeeze(pocFluxMonthlyDhNbyMethod(:,5,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;
fracHausgartenDataPointsDhAnalysis = (sum(squeeze(pocFluxMonthlyDhNbyMethod(:,6,:,:)),'all','omitnan')/nObsForDhAnalysis)*100;

fprintf('\n%0.1f%% data points at BATS/OFP,', fracBatsDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at HOT/ALOHA, and', fracHotalohaDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at OSP,', fracOspDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at PAP-SO,', fracPapsoDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at HAUSGARTEN.', fracHausgartenDataPointsDhAnalysis)
fprintf('\n%0.1f%% data points at EqPac,', fracEqpacDataPointsDhAnalysis)

nRadionuclideDataPoints = sum(squeeze(pocFluxMonthlyDhNbyMethod(:,:,:,2)),'all','omitnan');
fracRadionuclideDataPoints = (nRadionuclideDataPoints/nObsForDhAnalysis)*100;
fprintf('\nWe have %0.1f%% radionuclide data points for summaries by depth horizon.', fracRadionuclideDataPoints)

% Compare with length of the raw data set
nEntriesRawDatasetByLoc = zeros(NUM_LOCS,1);
for iLoc = 1:NUM_LOCS
    currStatData = D.POC_mmol_m2_d(D.tag == STATION_NAMES{iLoc},:);
    nEntriesRawDatasetByLoc(iLoc) = nnz(~isnan(currStatData));
end
nEntriesRawDatasetTotal = sum(nEntriesRawDatasetByLoc); % = nnz(~isnan(dataRaw.POC_mmol_m2_d));
fprintf('\nWe have %d data points IN TOTAL, of which', nEntriesRawDatasetTotal)

fracBatsDataPointsRaw   = (nEntriesRawDatasetByLoc(4)/nEntriesRawDatasetTotal)*100;
fracOspDataPointsRaw    = (nEntriesRawDatasetByLoc(2)/nEntriesRawDatasetTotal)*100;
fracPapsoDataPointsRaw  = (nEntriesRawDatasetByLoc(3)/nEntriesRawDatasetTotal)*100;
fracEqpacDataPointsRaw  = (nEntriesRawDatasetByLoc(1)/nEntriesRawDatasetTotal)*100;
fracHotalohaDataPointsRaw  = (nEntriesRawDatasetByLoc(5)/nEntriesRawDatasetTotal)*100;
fracHusgartenDataPointsRaw  = (nEntriesRawDatasetByLoc(6)/nEntriesRawDatasetTotal)*100;

fprintf('\n%0.1f%% data points at BATS/OFP', fracBatsDataPointsRaw)
fprintf('\n%0.1f%% data points at HOT/ALOHA', fracHotalohaDataPointsRaw)
fprintf('\n%0.1f%% data points at OSP', fracOspDataPointsRaw)
fprintf('\n%0.1f%% data points at PAP-SO', fracPapsoDataPointsRaw)
fprintf('\n%0.1f%% data points at HAUSGARTEN', fracHusgartenDataPointsRaw)
fprintf('\n%0.1f%% data points at EqPac', fracEqpacDataPointsRaw)

nEntriesRawDatasetRadionuclideTotal = length(D.POC_mmol_m2_d(strcmp(D.method,'radionuclide'),:));
fracRadionuclidesDatasetRaw = (nEntriesRawDatasetRadionuclideTotal/nEntriesRawDatasetTotal)*100;
fprintf('\nWe have %0.1f%% radionuclide data points in total.\n', fracRadionuclidesDatasetRaw)

% Print number of data points for reference
tablendp = zeros(NUM_TARGET_DEPTHS,NUM_LOCS,2);
for iLoc = 1:NUM_LOCS
    for iDh = 1:NUM_TARGET_DEPTHS
        for iMethod = 1:2
            ndp = sum(squeeze(pocFluxMonthlyDhNbyMethod(:,iLoc,iDh,iMethod)),'omitnan');
            tablendp(iDh,iLoc,iMethod) = ndp;
        end
    end
end
checkNumEntries = tablendp(1,6,2); % see in 6th loc, zeu for method 2 (radionuclide)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 9 - SAVE THE DATA
% -------------------------------------------------------------------------

save(fullfile('.','data','processed',filenameOutputFluxCompilation),...
    'pocFluxRawProfileValues_cell','pocFluxRawProfileErrRand_cell',...
    'pocFluxRawProfileErrSys_cell','pocFluxRawProfileDepths_cell',...
    'pocFluxRawProfileDataType_cell','pocFluxRawProfileValues',...
    'pocFluxRawProfileErrRand','pocFluxRawProfileErrSys',...
    'pocFluxRawProfileDepths','pocFluxRawProfileDataType',...
    'pocFluxRawDhValues_cell','pocFluxRawDhDepths_cell',...
    'pocFluxRawDhTag_cell','pocFluxRawDhDataType_cell',...
    'pocFluxMonthlyProfileAvg','pocFluxMonthlyProfileErrTot',...
    'pocFluxMonthlyProfileN','pocFluxMonthlyProfileDepths',...
    'pocFluxMonthlyDhAvg','pocFluxMonthlyDhErrTot','pocFluxMonthlyDhN',...
    'pocFluxAnnualDhAvg','pocFluxAnnualDhErrTot','pocFluxAnnualDhN',...
    'pocFluxAnnualProfileAvg','pocFluxAnnualProfileErrTot',...
    'pocFluxAnnualProfileN','pocFluxAnnualProfileDepths',...
    'TRAPRAD_TABLE_POC','LOC_DEPTH_HORIZONS_POC_OBS_ADAPT',...
    'LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED','-v7.3')

fprintf('\nThe POC flux compilation data have been saved correctly.\n')