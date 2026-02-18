function createSlamsInputData(fullpathRawDataDir, fullpathModelInputDataDir,...
   fullpathInterimDataDir, filenameInputPar0, filenameInputNpp, filenameInputChla,...
   filenameInputRho, filenameInputOmegaCalc, filenameInputNit, filenameInputSil,...
   filenameInputPhos, filenameInputOxy, filenameInputTemp, filenameInputDust,...
   filenameInputMld, filenameInputMesozoo, filenameInputDynVisco,...
   filenameInputNumDaylightHours, finenameInputMask, filenameInputRunGrid,...
   filenameInputNumDepthLayers, choiceOutputToForceModel, choiceLightDistribution,...
   choiceTypeGridDomain)

% CREATESLAMSINPUTDATA Pre-processing tool for generating SLAMS input data files.
% This script prepares input data from external data repositories based on 
% a specified grid domain, fills NaN values strategically (which are critical 
% to handle as SLAMS cannot run with NaNs), formats the output and organises 
% it into run instance folders required for simulation execution.
%
%   Version 1.0 - Completed 10 Jan 2026   
%
% Folder paths:
%   fullpathRawDataDir        = './data/raw/';
%   fullpathModelInputDataDir = './tests/LOCALTS6/modelinputdata/';
% 
% Input:
%   filenameInputPar0             = 'par0_modis.mat';
%   filenameInputNpp              = 'npp_bicep.mat';
%   filenameInputChla             = 'chla_occci.mat';
%   filenameInputRho              = 'rho_calculated_woa23.mat';
%   filenameInputOmegaCalc        = 'omegacalcite_co2sys.mat';
%   filenameInputNit              = 'nit_monthly_woa23.mat';
%   filenameInputSil              = 'sil_monthly_woa23.mat';
%   filenameInputPhos             = 'phos_monthly_woa23.mat';
%   filenameInputOxy              = 'oxy_monthly_woa23.mat';
%   filenameInputTemp             = 'temp_monthly_woa23.mat';
%   filenameInputDust             = 'dustflux_cmip6_ncarcesm2.mat';
%   filenameInputMld              = 'mld_ifremer.mat';
%   filenameInputMesozoo          = 'mesozoo_cmip6_pisces.mat';
%   filenameInputDynVisco         = 'dynamicvisco_calculated_woa23.mat';
%   filenameInputNumDaylightHours = 'ndaylighthours_daily_calculated.mat';
%   finenameInputMask             = 'mask_custom_icefrac_cmems_chla_occci.mat';
%   filenameInputRunGrid          = 'grid_run.mat';
%   filenameInputNumDepthLayers   = 'waterColNumDepthLayers.txt';
% 
% Choices:
%     choiceOutputToForceModel = false; 
%       true : outputs are .bin files to force the model
%       false: outputs are .mat files for visualisation
%     choiceLightDistribution = 1; 
%       1: light according to photoperiod distribution (1 timestep minimum)
%       2: light in two time steps (out of three)
%       3: always light
%     choiceTypeGridDomain = 1; 
%       1: global
%       2: local
%       3: global with zoom in factor
% 
% addpath(fullpathRawDataDir);
% addpath(fullpathModelInputDataDir);
% addpath(genpath('./modelresources/'));
% addpath(genpath('./code/matlab/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Load configuration parameters used in the model
config = loadModelConfigurationParameters(fullpathModelInputDataDir,...
    filenameInputRunGrid,filenameInputNumDepthLayers,choiceTypeGridDomain);

% Output filenames
filenameOutputLog          = 'logCreateInputData.txt';     % always an output
filenameOutputLonlatRef    = 'lonlatRef.txt';              % always an output
filenameOutputForcingArray = 'slamsForcing.mat';           

% Datasets metadata, all data already comes lat x lon x (depth) x time
datasetsMetadata = { 
  % {dataset name,  input filename,        variable names,                                             output filename,                  hasDepth?, suffersFromIceCover?, needsDensityAdjustment?, needsUnitsConversion?, unitFactor}
    'Par0',     filenameInputPar0,        {'par0_lat','par0_lon','par0'},                              config.filenameModelPar0,         false, true,  false, false, NaN; % W m-2
    'Npp',      filenameInputNpp,         {'npp_lat','npp_lon','npp_avg'},                             config.filenameModelNpp,          false, true,  false, true,  1/(config.molarMassCarbon*1e3); % mg C m-2 d-1 --> mol C m-2 s-1 (another units conversion is applied later to convert from d-1 to s-1)
    'Chla',     filenameInputChla,        {'chla_lat','chla_lon','chla'},                              config.filenameModelChla,         false, true,  false, false, NaN; % mg chla m-3
    'Rho',      filenameInputRho,         {'rho_lat','rho_lon','rho_depth','rho'},                     config.filenameModelRho,          true,  false, false, true,  1/1e3; % kg m-3 --> g cm-3
    'OmegaC',   filenameInputOmegaCalc,   {'glodap_lat','glodap_lon','glodap_depth','omegacalcite'},   config.filenameModelOmegaCalc,    true,  false, false, false, NaN; % unitless
    'Nit',      filenameInputNit,         {'woa_lat','woa_lon','woa_depth_nit','nit'},                 config.filenameModelNit,          true,  false, true,  false, NaN; % umol kg-1 --> umol L-1 (=mmol m-3) 
    'Sil',      filenameInputSil,         {'woa_lat','woa_lon','woa_depth_sil','sil'},                 config.filenameModelSil,          true,  false, true,  false, NaN; % umol kg-1 --> umol L-1 (=mmol m-3) 
    'Phos',     filenameInputPhos,        {'woa_lat','woa_lon','woa_depth_phos','phos'},               config.filenameModelPhos,         true,  false, true,  false, NaN; % umol kg-1 --> umol L-1 (=mmol m-3) 
    'Oxy',      filenameInputOxy,         {'woa_lat','woa_lon','woa_depth_oxy','oxy'},                 config.filenameModelOxy,          true,  false, true,  true,  config.molarVolumeOxygen/1e3; % umol kg-1 --> mL L-1
    'Temp',     filenameInputTemp,        {'woa_lat','woa_lon','woa_depth_temp','temp'},               config.filenameModelTemp,         true,  false, false, false, NaN; % ºC
    'Dust',     filenameInputDust,        {'dustflux_lat','dustflux_lon','dustflux'},                  config.filenameModelDust,         false, true,  false, true,  config.dustToClayRatio/(24*3600*1e3); % mg dust m-2 d-1 --> g clay m-2 s-1
    'Mld',      filenameInputMld,         {'mld_lat','mld_lon','mld'},                                 config.filenameModelMld,          false, false, false, false, NaN; % m
    'Mesozoo',  filenameInputMesozoo,     {'mesozoo_lat','mesozoo_lon','mesozoo_depth','mesozoo'},     config.filenameModelMesozoo,      true,  false, false, false, NaN; % mg C m-3
    'DynVisco', filenameInputDynVisco,    {'dynvisco_lat','dynvisco_lon','dynvisco_depth','dynvisco'}, config.filenameModelDynVisco,     true,  false, false, true,  1e3/1e2; % kg m-1 s-1 --> g cm-1 s-1
};

accessoryDatasetsMetadata = {
    'Ndh',     filenameInputNumDaylightHours, {'par0_lat','par0_lon','par0daylighthours'},           [],                              false, false, false, false, NaN;
};

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - TIME STEP CONFIGURATION
% -------------------------------------------------------------------------

if choiceOutputToForceModel
    nCalendarTimeSteps = 365; % daily interpolation
    nModelledTimeSteps = nCalendarTimeSteps*config.nTimestepsPerDay;
else
    nCalendarTimeSteps = 12; % monthly interpolation
    nModelledTimeSteps = nCalendarTimeSteps;
end

qTime = linspace(1,12,nCalendarTimeSteps)'; % month fraction    
timeStepLengthInHours = 24 / config.nTimestepsPerDay;

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - OPEN OUTPUT FILES
% -------------------------------------------------------------------------

% Lon/lat positioning
lonlatRefID = fopen(fullfile(fullpathModelInputDataDir,filenameOutputLonlatRef),'w');

% Depth layer count
numdepthlayersID = fopen(fullfile(fullpathModelInputDataDir,filenameInputNumDepthLayers),'w'); 

% Log progress
logID = fopen(fullfile(fullpathModelInputDataDir,filenameOutputLog),'w'); 

% Write to log
fprintf(logID, 'Loop through %d potential locations \n', config.nLocs);
fprintf(logID, 'Lon. increment: %.6f deg, Lat. increment: %.6f deg\n', config.lonIncrement, config.latIncrement);
if choiceTypeGridDomain == 1 || choiceTypeGridDomain == 3 % global
    if config.isLonUniform && config.isLatUniform
        fprintf(logID, 'The grid spacing is uniform.\n');
    else
        fprintf(logID, 'Warning: The grid spacing is not uniform. Stopping execution.\n');
        fclose(logID); % ensure the log file is closed before exiting
        error('Grid spacing is not uniform. Execution terminated.');
    end
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - DEFINE INTERPOLANTS
% -------------------------------------------------------------------------

% Structure to store interpolants to regrid oceanographic data from
% default dataset grid to chosen model grid
interpolantStruct = struct(); 
interpolantStruct = createDatasetInterpolant(interpolantStruct,datasetsMetadata,...
    fullpathRawDataDir,finenameInputMask,filenameInputRho,logID); % main metadata
interpolantStruct = createDatasetInterpolant(interpolantStruct,accessoryDatasetsMetadata,...
    fullpathRawDataDir,finenameInputMask,filenameInputRho,logID); % accessory metadata

% Save in case we don't want to run the whole process again
% save(fullfile('data','interim','inputDataInterpolants.mat'),'interpolantStruct','-v7.3')
% load(fullfile(fullpathInterimDataDir,'inputDataInterpolants.mat'),'interpolantStruct')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - DEFINE ARRAY TO VISUALISE OUTPUTS
% -------------------------------------------------------------------------

% Structure to store outputs for visualisation
mapStruct = initialiseMapStruct(datasetsMetadata,config.maxNoDepthLayers,...
    nModelledTimeSteps,config.nLocs);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - EXTRACT OCEANOGRAPHIC DATA FOR EACH CELL IN THE GLOBAL GRID
% -------------------------------------------------------------------------

% Minimising the number of tasks within the main loop over simulation locations
% improves performance. Many tasks have therefore been moved to the
% pre-processing function `createDatasetInterpolant` to optimise efficiency.

tic; % start timing the processing

for iLoc = 1:config.nLocs 
    fprintf(logID, '%d\n', iLoc);

    % Extract location-specific data
    [qLat,qLon,qDepth,nDepthLayers] = extractLocationData(config.gridLats,...
        config.gridLons,config.gridDepths,iLoc,logID);
    
    % Calculate depth boundaries
    [zUppBounds,zLowBounds] = calculateDepthBoundaries(qDepth,config.defaultDepthLayerThickness);

    % Initialise structure to store processed data for the current location
    processedData = struct();
    
    for iDataset = 1:height(datasetsMetadata)
        % Extract dataset-specific information required for interpolation
        [datasetName,hasDepth,F,lat,lon] = prepareInterpolation(...
            datasetsMetadata,iDataset,interpolantStruct,fullpathRawDataDir);
        %fprintf(logID, '%s\n', datasetName);
 
        % Perform interpolation of the dataset to the current location's coordinates
        data = interpolateDataset(F,qLat,qLon,qDepth,qTime,hasDepth);
  
        % Handle missing (NaN) values, excluding datasets where NaNs have 
        % been already taken care of, suffersFromIceCover = true)
        if ~ismember(datasetName,{'Par0','Npp','Chla','Dust'}) 
            data = manageGapsInOceanData(data,qDepth,qTime,hasDepth,F,qLat,qLon,...
                config.lonIncrement,config.latIncrement,lat(end),...
                config.maxIterForNeighbourhoodExpansion,logID);

            % Ensure no negative values (except for temperature datasets)
            if ~strcmp(datasetName,'Temp') && any(data(:) < 0) 
                throwError(logID,'Negative values detected in %s dataset \n',datasetName);
            end
            % Verify that there are no unresolved NaN values
            if any(isnan(data(:)))
                throwError(logID,'NaN values detected in %s dataset \n',datasetName);
            end
            
        % Handle datasets affected by ice cover, where NaNs result from masking
        % (periods without photosynthesis or aeolian dust flux deposition). 
        % We cannot pass the Fortran code forcing data with NaNs so we should 
        % replace those with 0    
        else   
            if any(data(:) < 0) 
                throwError(logID,'Negative values detected in %s dataset.\n',datasetName);
            end
            if choiceOutputToForceModel && any(isnan(data(:)))
                fprintf(logID,'NaN values detected in %s dataset, replacing with 0.\n',datasetName);
                data(isnan(data)) = 0; 
            end
        end

        % Shape output if required
        if ~ismember(datasetName, {'Npp','Par0'}) % exclude specific variables
            data = shapeOutputArray(data,choiceOutputToForceModel,hasDepth,config.nTimestepsPerDay);
        end
        
        % Store processed data
        processedData.(datasetName) = data;

    end % iDataset
 
    % Resolve interdependencies between datasets (NPP and PAR0)
    [processedData,avgPar0daylight] = processLightAndNPP(processedData,...
        interpolantStruct,qLat,qLon,accessoryDatasetsMetadata,choiceLightDistribution,...
        timeStepLengthInHours,config.nTimestepsPerDay,choiceOutputToForceModel);

    % Save outputs
    mapStruct = saveOutputs(mapStruct,iLoc,qLat,qLon,nDepthLayers,processedData,...
        choiceOutputToForceModel,fullpathModelInputDataDir,lonlatRefID,numdepthlayersID,...
        zUppBounds,zLowBounds,config.filenameModelZub,config.filenameModelZlb,...
        datasetsMetadata,avgPar0daylight);
    
    clear data processedData
  
end % iLoc
    
% Finalise and close files
% mapStruct IS GOING TO CONTAIN NANS IF THIS SCRIPT IS RUNNING TO CREATE
% .BIN FILES AS OPOSED TO VISUALISE INPUT
save(fullfile(fullpathModelInputDataDir,filenameOutputForcingArray),'mapStruct','-v7.3') 
fprintf(logID, 'Processing completed in %.2f seconds.\n', toc);    
fclose(logID);
fclose(lonlatRefID);
fclose(numdepthlayersID);

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function interpolantStruct = createDatasetInterpolant(interpolantStruct,...
    datasetsMetadata,fullpathInputDataDir,finenameInputMask,filenameInputRho,logID)

    for iDataset = 1:size(datasetsMetadata,1)

        % Extract dataset metadata
        [datasetName,filename,listFilenameVars,outputFilename,...
            hasDepth,suffersFromIceCover,needsDensityAdjustment,...
            needsUnitsConversion,unitFactor] = datasetsMetadata{iDataset,:};
    
        % Load dataset
        dataStruct = load(fullfile(fullpathInputDataDir,filename),listFilenameVars{:});

        % Map dataset variables dynamically (assuming consistent naming conventions)
        lat = dataStruct.(listFilenameVars{1});
        lon = dataStruct.(listFilenameVars{2});
        if hasDepth
            depth = dataStruct.(listFilenameVars{3});
            data = dataStruct.(listFilenameVars{4});
        else
            depth = [];
            data = dataStruct.(listFilenameVars{3});
        end

        % Handle gap-filling and masking for surface ocean data affected by 
        % ice cover
        if suffersFromIceCover
            if ~exist('mask', 'var') % load mask
                load(fullfile('data','raw',finenameInputMask),'mask','mask_lat','mask_lon')
                mask = single(mask);
            end
            data = manageGapsInOceanDataIceMasked(data,mask,lat,lon,mask_lat,mask_lon,(1:12),logID);
            clear mask
        end
        
        % Density adjustment
        if needsDensityAdjustment
            if ~exist('rho', 'var') % load density
                load(fullfile('data','raw',filenameInputRho),'rho','rho_lat','rho_lon','rho_depth')
            end
            data = applyDensityConversion(data,rho,lat,lon,depth,rho_lat,rho_lon,rho_depth,(1:12));
            clear rho
        end

        % Units conversion
        if needsUnitsConversion
            data = data .* unitFactor;
        end

        % Create interpolant for daily data (needs 365 points in time)
        if strcmp(datasetName,'Ndh')

            [X, Y, T] = ndgrid(lat, lon, (1:365)');
            interpolantStruct.(datasetName) = griddedInterpolant(X, Y, T, data, 'linear', 'linear');

            fprintf(logID, 'Interpolant for %s created \n',datasetName);
            clear X Y T data
                
        % Create interpolant for data that comes on a monthly format
        else

            if isempty(depth) % surface variables

                if ndims(data) == 2 % lat x lon
                    data = repmat(data, 1, 1, 12); % repeat across 12 months for surface data
                end
                [X, Y, T] = ndgrid(lat, lon, (1:12)');
                interpolantStruct.(datasetName) = griddedInterpolant(X, Y, T, data, 'linear', 'linear');

                fprintf(logID, 'Interpolant for 2D %s created \n',datasetName);
                clear X Y T data
                
            else % depth-dependent variables

                if ndims(data) == 3 % lat x lon x depth
                    data = repmat(data, 1, 1, 1, 12); % repeat across 12 months for depth data
                end
                [X, Y, Z, T] = ndgrid(lat, lon, depth, (1:12)');
                interpolantStruct.(datasetName) = griddedInterpolant(X, Y, Z, T, data, 'linear', 'linear');

                fprintf(logID, 'Interpolant for 3D %s created \n',datasetName);
                clear X Y T Z data
                
            end
        end

    end % iDataset
    
    clear dataStruct
    
end % createDatasetInterpolant
    
% *************************************************************************

function mapStruct = initialiseMapStruct(datasetsMetadata,maxNumDepthLayers,...
    nModelledTimeSteps,nLocs)

    mapStruct = struct(); % initialise
    
    % Data with monthly output
    for i = 1:size(datasetsMetadata,1)
        datasetName = datasetsMetadata{i,1};
        hasDepth = datasetsMetadata{i,5};
        if ~strcmp(datasetName,'Par0') % exclude specific variables
            if hasDepth
                mapStruct.(['monthly' datasetName]) = NaN(maxNumDepthLayers,nModelledTimeSteps,nLocs);
            else
                mapStruct.(['monthly' datasetName]) = NaN(nModelledTimeSteps,nLocs);
            end
        end
    end
    
    % Data with daily output
    mapStruct.dailyPar0         = NaN(365,nLocs); % for plotting
    mapStruct.dailyPar0daylight = NaN(365,nLocs); % to calculate phytoplankton probabilities
   
end % initialiseMapStruct

% *************************************************************************

function [qLat,qLon,qDepth,nDepthLayers] = extractLocationData(gridLats,...
    gridLons,gridDepths,iLoc,logID)

    qLat = gridLats(iLoc);
    qLon = gridLons(iLoc);
    nDepthLayers = sum(~isnan(gridDepths(:,iLoc))); % count valid depth layers
    qDepth = gridDepths(1:nDepthLayers,iLoc); % extract valid depths
    
    % Log progress
    logEntry = sprintf('[%s] Location %d: lat/lon %.2f/%.2f, depthLayers %d\n',...
        datestr(now), iLoc, qLat, qLon, nDepthLayers);
    fprintf(logID,logEntry);
    
end % extractLocationData

% *************************************************************************

function [upperBounds,lowerBounds] = calculateDepthBoundaries(depths,maxThickness)
    
    halfThickness = maxThickness / 2;
    upperBounds = depths - halfThickness;
    lowerBounds = depths + halfThickness;
    
end % calculateDepthBoundaries

% *************************************************************************

function [datasetName,hasDepth,F,lat,lon] = prepareInterpolation(...
    datasetsMetadata,iDataset,interpolantStruct,fullpathInputDataDir)

    % Extract dataset metadata
    [datasetName,filename,listFilenameVars,outputFilename,...
        hasDepth,suffersFromIceCover,needsDensityAdjustment,...
        needsUnitsConversion,unitFactor] = datasetsMetadata{iDataset,:};

    % Access dataset interpolant from struct
    F = interpolantStruct.(datasetName); 

    % Map variables dynamically (assuming consistent naming conventions)
    dataStruct = load(fullfile(fullpathInputDataDir,filename),listFilenameVars{1:2});
    lat = dataStruct.(listFilenameVars{1});
    lon = dataStruct.(listFilenameVars{2});
        
end % prepareInterpolation

% *************************************************************************

function interpolatedData = interpolateDataset(F,qLat,qLon,qDepth,qTime,hasDepth)

    % Query grid construction
    [qX, qY, qT] = ndgrid(qLat,qLon,qTime'); % query grid
    if hasDepth
        [qX, qY, qZ, qT] = ndgrid(qLat,qLon,qDepth',qTime');
        interpolatedData = squeeze(F(qX, qY, qZ, qT)); % query interpolant with the grid
    else
        interpolatedData = squeeze(F(qX, qY, qT));
    end

end % interpolateDataset

% *************************************************************************

function outputArray = shapeOutputArray(inputArray,choiceOutputToForceModel,...
    hasDepth,nTimestepsPerModelledDay)
    
    if choiceOutputToForceModel
        % Repeat array for forcing model with timestep resolution
        if hasDepth % 2D array
            outputArray = repelem(inputArray,1,nTimestepsPerModelledDay); 
        else % 1D case
            outputArray = repelem(inputArray,nTimestepsPerModelledDay); 
        end
    else
        % Output monthly data
        outputArray = inputArray; 
    end
    
end % shapeOutputArray

% *************************************************************************

function mapStruct = saveOutputs(mapStruct,iLoc,qLat,qLon,nDepthLayers,...
    processedData,choiceOutputToForceModel,fullpathOutputDataDir,lonlatRefID,...
    numdepthlayersID,zUppBounds,zLowBounds,filenameOutputZub,filenameOutputZlb,...
    datasetsMetadata,avgPar0daylight)
    
    % Print to output file
    fprintf(lonlatRefID, '%d, %.2f, %.2f\n', iLoc, qLat, qLon);

    if choiceOutputToForceModel % save binary output for forcing model

        % Print to output file
        fprintf(numdepthlayersID, '%d\n', nDepthLayers);
        
        % Create output directory
        locOutputDir = fullfile(fullpathOutputDataDir, sprintf('run_%d/', iLoc));
        mkdir(locOutputDir);

        % Forcing data
        for i = 1:size(datasetsMetadata,1)
            datasetName = datasetsMetadata{i,1};
            outputFilename = datasetsMetadata{i,4};
            outputVar = processedData.(datasetName);

            % Save to binary files 
            fid = fopen(fullfile(locOutputDir,outputFilename),'w');
            if fid < 0
                error('Failed to open file for writing: %s', outputFilename);
            end
            fwrite(fid, outputVar, 'float64', 'ieee-be');
            fclose(fid);
            
            % ACTIVATE IF RUNNING LOCALLY AND WISH TO SEE INPUT DATA
            %% Save to matlab structure
            %hasDepth = datasetsMetadata{i,5};
            %if ~strcmp(datasetName, 'Par0')
            %    fieldName = ['monthly' datasetName];
            %    if hasDepth
            %       mapStruct.(fieldName)(1:nDepthLayers,:,iLoc) = processedData.(datasetName);
            %    else
            %        mapStruct.(fieldName)(:,iLoc) = processedData.(datasetName);
            %    end
            %end
        end

		% ACTIVATE IF RUNNING LOCALLY AND WISH TO SEE INPUT DATA
        % Additional forcing data assignment to matlab structure
        %mapStruct.dailyPar0daylight(:,iLoc) = avgPar0daylight; % to calculate phytoplankton probabilities

        % Binary files for depth layer data
        fid = fopen(fullfile(locOutputDir,filenameOutputZub),'w');
        if fid < 0
            error('Failed to open file for writing: %s', filenameOutputZub);
        end
        fwrite(fid, zUppBounds, 'float64', 'ieee-be');
        fclose(fid);

        fid = fopen(fullfile(locOutputDir, filenameOutputZlb), 'w');
        if fid < 0
            error('Failed to open file for writing: %s', filenameOutputZlb);
        end
        fwrite(fid, zLowBounds, 'float64', 'ieee-be');
        fclose(fid);
        
    else % update mapStruct
 
        % Forcing data
        for i = 1:size(datasetsMetadata,1)
            datasetName = datasetsMetadata{i,1};
            hasDepth = datasetsMetadata{i,5};

            % Exclude specific variables by name
            if ~strcmp(datasetName, 'Par0')
                % Handle datasets with or without depth
                fieldName = ['monthly' datasetName];
                if hasDepth
                    mapStruct.(fieldName)(1:nDepthLayers,:,iLoc) = processedData.(datasetName);
                else
                    mapStruct.(fieldName)(:,iLoc) = processedData.(datasetName);
                end
            else
                % Special case for 'Par0' variable
                mapStruct.(['daily', datasetName])(:,iLoc) = processedData.(datasetName); % =avgPar0daylightTimeStepAdapted (for plotting)
            end
        end

        % Additional data assignment for avgPar0daylight
        mapStruct.dailyPar0daylight(:,iLoc) = avgPar0daylight; % to calculate phytoplankton probabilities

    end
    
end % saveOutputs

% *************************************************************************

function [processedData,avgPar0daylight] = processLightAndNPP(processedData,...
    interpolantStruct,qLat,qLon,accessoryDatasetsMetadata,choiceLightDistribution,...
    timeStepLengthInHours,nTimeStepsModelledDay,choiceOutputToForceModel)
   
    npp = processedData.('Npp');
    par0 = processedData.('Par0');

    % Interpolate number of daylight hours at qLat and qLon
    datasetName = accessoryDatasetsMetadata{1}; % extract metadata for number of daylight hours
    F = interpolantStruct.(datasetName);
    [qX, qY, qT] = ndgrid(qLat,qLon,(1:365)'); % query grid
    qNumDaylightHours = squeeze(F(qX, qY, qT));

    % Daily average PAR0 interpolation or direct assignment
    if ~choiceOutputToForceModel
        % From monthly to daily resolution
        midMonthDay = [15, 46, 77, 106, 137, 167, 198, 228, 259, 289, 320, 350]'; % define the middle day of each month (15th day of each month)
        F = griddedInterpolant(midMonthDay, par0, 'linear', 'linear');
        avgPar0daylight = F((1:365)'); % W m-2
    else
        % All good
        avgPar0daylight = par0; % W m-2
    end
    
    % Sum over number of daylight hours
    totPar0daylight = avgPar0daylight.*qNumDaylightHours*3600; % J m-2

    % Determine time steps with light based on user choice
    switch choiceLightDistribution
        case 1
            nTimeStepsWithLight = max(round(qNumDaylightHours./timeStepLengthInHours),1);
        case 2 
            nTimeStepsWithLight = 2 * ones(365,1); % only two 
        case 3
            nTimeStepsWithLight = 3 * ones(365,1); % three      
    end
    nSunnySecondsPerDay = (nTimeStepsWithLight./nTimeStepsModelledDay).*3600.*24;

    % Compute average PAR0 per time step (adapted for daylight hours)
    avgPar0daylightTimeStepAdap = totPar0daylight ./ ...
        (nTimeStepsWithLight * timeStepLengthInHours * 3600); % W m-2

    % Prepare output arrays for PAR0
    if choiceOutputToForceModel
        forcPar0 = zeros(365,nTimeStepsModelledDay); % W m-2
        for iDay = 1:365
            forcPar0(iDay,1:nTimeStepsWithLight(iDay)) = ...
                repelem(avgPar0daylightTimeStepAdap(iDay),1,nTimeStepsWithLight(iDay));
        end
        forcPar0 = reshape(forcPar0', (365 * nTimeStepsModelledDay), 1);
    else
        forcPar0 = avgPar0daylightTimeStepAdap;
    end

    % Process NPP
    if choiceOutputToForceModel
        npp = npp./nSunnySecondsPerDay; % mol C m-2 d-1 --> mol C m-2 s-1
        forcNpp = zeros(365,nTimeStepsModelledDay); % mol C m-2 s-1
        for iDay = 1:365
            forcNpp(iDay,1:nTimeStepsWithLight(iDay)) = ...
                repelem(npp(iDay),1,nTimeStepsWithLight(iDay));
        end
        forcNpp = reshape(forcNpp', (365 * nTimeStepsModelledDay), 1);
    else
        npp = npp./(24*3600); % mol C m-2 d-1 --> mol C m-2 s-1        
        forcNpp = npp; 
    end
    
    % Save newly calculated NPP and PAR0 to processedData structure
    processedData.('Par0') = forcPar0; % time step adapted
    processedData.('Npp') = forcNpp; % time step adapted
 
end % processLightAndNPP

% *************************************************************************
        
% function qNpp = interpolateNppFromInsituObservations(qLat,qLon,fullpathDataDir,...
%     filenameNppObservations)
% 
%     if (qLat == 50 && qLon == -145) % OSP
% 
%         nppOSPobservationsID = fopen(strcat(fullpathDataDir,filenameNppObservations{1}),'r');    
%         nppObs = fscanf(nppOSPobservationsID,'%d\n'); 
%         fclose(nppOSPobservationsID);        
%         nppObs(nppObs < 0) = 0;
%         qNppObs = cleverTimeInterpolationLoc(nppObs,(1:12));
%         disp('Observed NPP for OSP:')
%         disp(qNppObs)
%         qNpp(:) = interp1((1:12), qNppObs, qTimes, 'linear'); % mg C m-2 d-1
% 
%     elseif (qLat == 31.6 && qLon == -64.2) % BATS
% 
%         nppBATSobservationsID = fopen(strcat(fullpathDataDir,filenameNppObservations{2}),'r');    
%         nppObs = fscanf(nppBATSobservationsID,'%d\n'); 
%         fclose(nppBATSobservationsID);
%         disp('Observed NPP for BATS/OFP:')
%         disp(nppObs)
%         qNpp(:) = interp1((1:12), nppObs, qTimes, 'linear'); % mg C m-2 d-1
% 
%     end
% 
% end % interpolateNppFromInsituObservations

end % createSlamsInputData