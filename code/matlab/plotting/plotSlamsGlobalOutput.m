
% ======================================================================= %
%                                                                         %
% This script extracts relevant SLAMS output data and processes it for    %
% later plotting. It requires the following scripts to have been run:     %
%   mapSlamsGlobalInput.m: produces slamsForcingGeo.mat                   %
%   readSlamsOutput.m: produces runsoutput.mat                            %
%   calculateBcpMetricsFromSlamsOutput.m: produces bcpmetrics.mat         %
%   processParticleNumberFromUvp.m: produces pocflux_bisson_16sc.mat      %
%                                                                         %
%   Version 1.0 - Completed 28 Feb 2025                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath(fullfile('.','modelresources','external')))
addpath(genpath(fullfile('.','modelresources','internal'))) 
addpath(genpath(fullfile('.','data','processed'))) 
addpath(genpath(fullfile('.','data','raw'))) 
addpath(genpath(fullfile('.','code','matlab'))) 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Name of experiment
testrunDir = 'GLOBALSTD';
figureSubfolder = 'GLOBALSTD';

% Choices to be applied
choiceTypeGridDomain = 1; % 1=global, 2=local, 3=global with zoom in factor

% Filename declarations
filenameOutputSlams           = 'runsoutput.mat';
filenameOutputSlamsBcpMetrics = 'bcpmetrics.mat';
filenameInputRunGrid          = 'grid_run.mat';
filenameInputNumDepthLayers   = 'waterColNumDepthLayers.txt';
filenameInputGeoForcing       = 'slamsForcingGeo.mat';
filenameInputNpp              = 'npp_bicep.mat';
filenameInputMask             = 'mask_custom_icefrac_cmems_chla_occci.mat';
filenameInputZeu              = 'zeu_calculated_chlaoccci_mldifremer_pointonepercentpar0.mat';

% Path declarations
fullpathTestDir           = fullfile('.','tests',testrunDir);
fullpathProcessedDataDir  = fullfile('.','data','processed');
fullpathInterimDataDir    = fullfile('.','data','interim');
fullpathModelInputDataDir = fullfile(fullpathTestDir,'modelinputdata');
fullpathModelRunsDir      = fullfile(fullpathTestDir,'modelruns');
fullpathConfigDir         = fullfile('.','config');

% Load configuration parameters used in the model
config = loadModelConfigurationParameters(fullpathModelInputDataDir,...
    filenameInputRunGrid,filenameInputNumDepthLayers,choiceTypeGridDomain);

% Log progress
logID = fopen(fullfile(fullpathModelRunsDir,'logCalculateBcpMetric.log'),'w');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD MODEL OUTPUT AND INPUT DATA
% -------------------------------------------------------------------------

% Load model output data
load(fullfile(fullpathModelRunsDir,filenameOutputSlams),'output') 

% % Load BCP metrics calculated
% load(fullfile(fullpathModelRunsDir,filenameOutputSlamsBcpMetrics),'metrics')

% Load model input data
load(fullfile(fullpathModelInputDataDir,filenameInputGeoForcing),'forcingDataGeo')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - PREPARE DATA FOR PLOTTING
% -------------------------------------------------------------------------

% Extract the output data that will be dealt with in this script (needs to
% be edited if more outputs needs to be plotted)
tmpOutput = extractModelOutput(output);

%%
% %%%% QUICK PLOTTING %%%%%
% geoPocFlux = NaN(config.nLats,config.nLons,config.maxNoSedTrapDeployDepths,12);
% 
% for iLoc = 1:config.nLocs
%     iLat = config.iyBb(iLoc);
%     iLon = config.ixBb(iLoc);
%     geoPocFlux(iLat,iLon,:,:) = squeeze(tmpOutput.monthlyOrgCarbonFlux(:,:,iLoc)); 
% end
% 
% myColourMap = config.mappingProps.myColourMapOceanVars;
% myData = geoPocFlux(:,:,36,6);
% lonVector = config.mappingProps.mapLons;
% latVector = config.mappingProps.mapLats;
% %%
% % Convert to gridded format
% [X, Y] = meshgrid(lonVector, latVector);
% 
% % Create figure and plot
% figure;
% pcolor(X, Y, myData);
% shading flat;  % Removes grid lines
% 
% % Customize color scaling and colormap
% caxis([0, 10]);  
% colormap(jet(100));
% colorbar;  % Add colorbar
% 
% % Set labels and title
% xlabel('Longitude');
% ylabel('Latitude');
% title('My Data Map');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - GET KEY DEPTH LEVELS FOR PLOTTING GLOBAL OUTPUT
% -------------------------------------------------------------------------

% Extract supporting data for calculations
[seqZeuMonthly,seqZeuAnnual,seqNppMonthly] = extractOutputSupportingData(...
    config,filenameInputZeu,filenameInputNpp,filenameInputMask,logID);

[keyDepthLevelsMonthlyData,keyDepthLevelsAnnualData] = ...
    defineKeyDepthsForPlottingGlobalOutput(config,seqZeuMonthly);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - COMPUTE UNCERTAINTIES FOR FLUXES AND PARTICLE NUMBERS
% -------------------------------------------------------------------------

[fluxSimulatedRelUncertaintyOverall,pNumSimulatedRelUncertaintyOverall] =...
    calculateBootstrapUncertainty(tmpOutput);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - EXTRACT DEPTH HORIZON DATA OF PARTICLE FLUXES, PARTICLE
% NUMBERS, PARTICLE ATTRIBUTES, SMS TERMS AND AUXILLIARY TERMS
% -------------------------------------------------------------------------

% Add fractional uncertainty to arrays and format them for further processing
targetOutput.flux.monthly.profile = prepareSeqData(tmpOutput,fluxSimulatedRelUncertaintyOverall,'flux');
targetOutput.pnum.monthly.profile = prepareSeqData(tmpOutput,pNumSimulatedRelUncertaintyOverall,'particleNumbers');
targetOutput.patt.annual.profile  = prepareSeqData(tmpOutput,[],'particleAttributes');
targetOutput.sms.annual.profile   = prepareSeqData(tmpOutput,[],'sms');
targetOutput.aux.annual.profile   = prepareSeqData(tmpOutput,[],'aux');

% Extract flux at key depth horizons
keyDepthLevelsMonthlyFluxData = repmat(keyDepthLevelsMonthlyData, 1, 1, 1, 1, 3); % repeat the 4x12x4448x2 array 3 times (3 tracers) along a new 5th dimension
[targetOutput.flux.monthly.dh,~,targetOutput.flux.annual.dh,~] =...
    extractKeyDepthHorizonDataInModelledFluxes(config,...
        targetOutput.flux.monthly.profile,keyDepthLevelsMonthlyFluxData);

% Compute flux material ratios at key depth horizons
[targetOutput.fluxrat.monthly.dh,targetOutput.fluxrat.annual.dh] =...
    calculateFluxRatios(config,targetOutput.flux.monthly.dh,targetOutput.flux.annual.dh);

% Extract particle numbers at key depth horizons and propagate uncertainty (in particle number units)
[targetOutput.pnum.monthly.dh,~,targetOutput.pnum.annual.dh,~] =...
    extractKeyDepthHorizonDataInModelledParticleNumbers(config,...
        targetOutput.pnum.monthly.profile,keyDepthLevelsMonthlyData);

% Extract particle attributes at key depth horizons
targetOutput.patt.annual.dh = extractKeyDepthHorizonDataInModelledParticleAttributes(...
    config,targetOutput.patt.annual.profile,keyDepthLevelsAnnualData);

% Extract SMS and aux terms at key depth horizons
[targetOutput.sms.annual.dh,targetOutput.aux.annual.dh] =...
    extractKeyDepthHorizonDataInModelledSmsAndAuxTerms(config,...
        targetOutput.sms.annual.profile,targetOutput.aux.annual.profile,...
        keyDepthLevelsAnnualData);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - COMPUTE EFFICIENCY METRICS (WITHOUT ERROR PROPAGATION)
% -------------------------------------------------------------------------

arrayFluxZeu   = reshape(targetOutput.flux.monthly.dh(1,:,:,1,1), 12, config.nLocs);
arrayFluxZmeso = reshape(targetOutput.flux.monthly.dh(3,:,:,1,1), 12, config.nLocs);
arrayNpp       = reshape(seqNppMonthly(:,:,1), 12, config.nLocs);

[teffAnnual,teffMonthly] = calculateEfficiencyMetricsWithoutErrorPropagation(...
    arrayFluxZmeso,arrayFluxZeu);

[peeffAnnual,peeffMonthly] = calculateEfficiencyMetricsWithoutErrorPropagation(...
    arrayFluxZeu,arrayNpp);

% Store results for Teff and PEeff using a struct update loop (allows
% adding Martin's b and z*)
effMetrics = {'teff','peeff'};
annualResults = {teffAnnual, peeffAnnual};
monthlyResults = {teffMonthly, peeffMonthly};
for i = 1:length(effMetrics)
    thisMetric = effMetrics{i};
    metrics.(thisMetric).monthly = monthlyResults{i}(:,:,:);
    metrics.(thisMetric).annual  = annualResults{i}(:,:);
end
                  
% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - PLOTS
% -------------------------------------------------------------------------

% plotBcpMetrics(figureSubfolder,'metric',metrics,config)
% plotFlux2Dslices(figureSubfolder,'flux',targetOutput,config)
plotAverageParticleProperties(figureSubfolder,'part',targetOutput,config)
% plotSmsAndAuxProperties(figureSubfolder,'suppl',targetOutput,config)
% plotContourPlotsForcingVsMetrics(figureSubfolder,'contour',metrics,forcingDataGeo,targetOutput,config)

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function [output,metrics,targetOutput] = extractTimeSeriesLocationsFromGlobalProcessedOutput(...
    config,output,metrics,targetOutput,fullpathConfigDir,filenameLocalLatitudes,...
    filenameLocalLongitudes)

    % Read the latitudes from the text file containing local data
    localLatitudes = readmatrix(fullfile(fullpathConfigDir,filenameLocalLatitudes));
    localLongitudes = readmatrix(fullfile(fullpathConfigDir,filenameLocalLongitudes));
    nLocalLocs = length(localLatitudes);

     % Initialize new indices array to extract locations
    selectedIndices = NaN(1,nLocalLocs);

    for iLoc = 1:nLocalLocs

        % Find the idx to the coordinates in the grid
        [~,iLon] = min(abs(config.lons(:) - localLongitudes(iLoc)));
        [~,iLat] = min(abs(config.lats(:) - localLatitudes(iLoc)));

        lonToSearch = config.lons(iLon);
        latToSearch = config.lats(iLat);

        isMatch = (config.gridLons == lonToSearch) & (config.gridLats == latToSearch);
        idxLoc = find(isMatch);

        if isempty(idxLoc)
            warning('No matching grid location found for Lat: %.2f, Lon: %.2f', latToSearch, lonToSearch);
            continue;
        end
        
        selectedIndices(iLoc) = idxLoc;
        
    end
    
    % Extract and overwrite original data using selected locations
    output.deployment.sedTrap = output.deployment.sedTrap(:,selectedIndices);
    output.deployment.imagSys = output.deployment.imagSys(:,selectedIndices);

    metrics.martinb.monthly = metrics.martinb.monthly(selectedIndices,:,:);
    metrics.martinb.annual = metrics.martinb.annual(selectedIndices,:);

    metrics.zstar.monthly = metrics.zstar.monthly(selectedIndices,:,:);
    metrics.zstar.annual = metrics.zstar.annual(selectedIndices,:);

    metrics.peeff.monthly = metrics.peeff.monthly(selectedIndices,:,:);
    metrics.peeff.annual = metrics.peeff.annual(selectedIndices,:);

    metrics.teff.monthly = metrics.teff.monthly(selectedIndices,:,:);
    metrics.teff.annual = metrics.teff.annual(selectedIndices,:);

    targetOutput.flux.monthly.profile = targetOutput.flux.monthly.profile(:,:,selectedIndices,:,:);
    targetOutput.pnum.monthly.profile = targetOutput.pnum.monthly.profile(:,:,:,selectedIndices,:);

end % extractTimeSeriesLocationsFromGlobalProcessedOutput

% *************************************************************************

function plotBcpMetrics(figureSubfolderName,figurePrefixName,metrics,config)

% Plot
fluxMetricsMonthlyMetadata = {...   
    'peeff',   squeeze(metrics.peeff.monthly(:,:)), 'Particle export efficiency',                   0, 0.70;
    'teff',    squeeze(metrics.teff.monthly(:,:)),  'Transfer efficiency (export depth to 1000 m)', 0, 0.30;
};

fluxMetricsAnnualMetadata = {...   
    'peeff',   squeeze(metrics.peeff.annual(:)), 'Particle export efficiency',                   0, 0.70;
    'teff',    squeeze(metrics.teff.annual(:)),  'Transfer efficiency (export depth to 1000 m)', 0, 0.30;
};

% Arrange into a structure the metadata
dataTypes = struct( ...
    'metadata', {fluxMetricsMonthlyMetadata, fluxMetricsAnnualMetadata}, ...
    'tag', {'output_global_monthly_flux_metrics', 'output_global_annual_flux_metrics'} ...
);

for iDatatype = 1:numel(dataTypes)
    datasetsMetadata = dataTypes(iDatatype).metadata;
    datasetsMetadata = unpackData(config,datasetsMetadata); % from location sequence structure to lat x lon array
    tagDataType = dataTypes(iDatatype).tag;
    
    % Initialise colour bar metadata
    caxisMin = zeros(height(datasetsMetadata),1);
    caxisMax = zeros(height(datasetsMetadata),1);
    for iDataset = 1:height(datasetsMetadata)
        caxisMin(iDataset) = datasetsMetadata{iDataset,4};
        caxisMax(iDataset) = datasetsMetadata{iDataset,5};
    end

    if iDatatype == 1 % monthly data (displayed as subplots in one figure)
        
        for iDataset = 1:height(datasetsMetadata)
            datasetName = lower(datasetsMetadata{iDataset,1});
            sgTitle = datasetsMetadata{iDataset,3};
            dataArray = datasetsMetadata{iDataset,2};
            plotOceanVariableMaps(dataArray,config.mappingProps.mapLons,...
                config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
                [],caxisMin(iDataset),caxisMax(iDataset),[],...
                config.mappingProps.labelMonths,...
                strcat(figureSubfolderName,'_',figurePrefixName,'_monthly_',tagDataType,'_',datasetName),...
                figureSubfolderName,sgTitle)
        end

    else % annual data (group into subplots)
        
        groupDatasets = [];
        titleStr = cell(height(datasetsMetadata),1);
        for iDataset = 1:height(datasetsMetadata)
            dataArray = datasetsMetadata{iDataset,2};
            groupDatasets = cat(3,groupDatasets,dataArray); % append
            titleStr{iDataset} = datasetsMetadata{iDataset,3};
        end

        isCommonColourBar = false;
        plotOceanVariableMaps(groupDatasets,config.mappingProps.mapLons,...
            config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
            [],caxisMin,caxisMax,isCommonColourBar,titleStr,...
            strcat(figureSubfolderName,'_',figurePrefixName,'_annual_',tagDataType),...
            figureSubfolderName,[])

    end
end % iDatatype

    % isAvgOrErr = 1; % 1 = average, 2 = upper boundary, 3 = lower boundary
    % tagStatTypes = {'avg', 'upper', 'lower'};
    % tagStatType = tagStatTypes{min(max(isAvgOrErr, 1), numel(tagStatTypes))};
    % 
    % fluxMetricsMonthlyMetadata = {...   
    %     'peeff',   squeeze(model.peeff.monthly(:,:,isAvgOrErr)),   'Particle export efficiency',                   0, 0.10;
    %     'teff',    squeeze(model.teff.monthly(:,:,isAvgOrErr)),    'Transfer efficiency (export depth to 1000 m)', 0, 0.10;
    %     'martinb', squeeze(model.martinb.monthly(:,:,isAvgOrErr)), 'Martin''s b',                                  0, 2;
    %     'zstar',   squeeze(model.zstar.monthly(:,:,isAvgOrErr)),   'Remineralisation length scale (m)',            100, 1500;
    % };
    % 
    % particleMetricsMonthlyMetadata = {...
    %     'xi_ex',   squeeze(model.xi.monthly(:,1,:,isAvgOrErr)), 'Slope PSD at export depth', -7, 7;
    %     'xi_mid',  squeeze(model.xi.monthly(:,2,:,isAvgOrErr)), 'Slope PSD at 500 m',        -7, 7;
    %     'xi_seq',  squeeze(model.xi.monthly(:,3,:,isAvgOrErr)), 'Slope PSD at 1000 m',       -7, 7;
    %     'xi_deep', squeeze(model.xi.monthly(:,4,:,isAvgOrErr)), 'Slope PSD at 2000 m',       -7, 7;
    % };
    % 
    % fluxMetricsAnnualMetadata = {...   
    %     'peeff',   squeeze(model.peeff.annual(:,isAvgOrErr)),   'Particle export efficiency',                   0, 0.10;
    %     'teff',    squeeze(model.teff.annual(:,isAvgOrErr)),    'Transfer efficiency (export depth to 1000 m)', 0, 0.10;
    %     'martinb', squeeze(model.martinb.annual(:,isAvgOrErr)), 'Martin''s b',                                  0, 2;
    %     'zstar',   squeeze(model.zstar.annual(:,isAvgOrErr)),   'Remineralisation length scale (m)',            100, 1500;
    % };
    % 
    % particleMetricsAnnualMetadata = {...
    %     'xi_ex',   squeeze(model.xi.annual(:,1,isAvgOrErr)), 'Slope PSD at export depth', -7, 7;
    %     'xi_mid',  squeeze(model.xi.annual(:,2,isAvgOrErr)), 'Slope PSD at 500 m',        -7, 7;
    %     'xi_seq',  squeeze(model.xi.annual(:,3,isAvgOrErr)), 'Slope PSD at 1000 m',       -7, 7;
    %     'xi_deep', squeeze(model.xi.annual(:,4,isAvgOrErr)), 'Slope PSD at 2000 m',       -7, 7;
    % };
    % 
    % % Arrange into a structure the metadata
    % dataTypes = struct( ...
    %     'metadata', {fluxMetricsMonthlyMetadata, particleMetricsMonthlyMetadata, fluxMetricsAnnualMetadata, particleMetricsAnnualMetadata}, ...
    %     'tag', {'output_global_monthly_flux_metrics', 'output_global_monthly_part_metrics', 'output_global_annual_flux_metrics', 'output_global_annual_part_metrics'} ...
    % );
    % 
    % for iDatatype = 1:numel(dataTypes)
    %     datasetsMetadata = dataTypes(iDatatype).metadata;
    %     datasetsMetadata = unpackData(config,datasetsMetadata); % from location sequence structure to lat x lon array
    % 
    %     tagDataType = dataTypes(iDatatype).tag;
    % 
    %     % Initialise colour bar metadata
    %     caxisMin = zeros(height(datasetsMetadata),1);
    %     caxisMax = zeros(height(datasetsMetadata),1);
    %     for iDataset = 1:height(datasetsMetadata)
    %         caxisMin(iDataset) = datasetsMetadata{iDataset,4};
    %         caxisMax(iDataset) = datasetsMetadata{iDataset,5};
    %     end
    % 
    %     if iDatatype <= 2 % monthly data (displayed as subplots in one figure)
    % 
    %         for iDataset = 1:height(datasetsMetadata)
    %             datasetName = lower(datasetsMetadata{iDataset,1});
    %             sgTitle = datasetsMetadata{iDataset,3};
    %             dataArray = datasetsMetadata{iDataset,2};
    %             plotOceanVariableMaps(dataArray,config.mappingProps.mapLons,...
    %                 config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
    %                 [],caxisMin(iDataset),caxisMax(iDataset),[],...
    %                 config.mappingProps.labelMonths,...
    %                 strcat(figurePrefixName,'_',tagDataType,'_',tagStatType,'_',datasetName),...
    %                 figureSubfolderName,sgTitle)
    %         end
    % 
    %     else % annual data (group into subplots)
    % 
    %         groupDatasets = [];
    %         titleStr = cell(height(datasetsMetadata),1);
    %         for iDataset = 1:height(datasetsMetadata)
    %             dataArray = datasetsMetadata{iDataset,2};
    %             groupDatasets = cat(3,groupDatasets,dataArray); % append
    %             titleStr{iDataset} = datasetsMetadata{iDataset,3};
    %         end
    % 
    %         isCommonColourBar = false;
    %         plotOceanVariableMaps(groupDatasets,config.mappingProps.mapLons,...
    %             config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
    %             [],caxisMin,caxisMax,isCommonColourBar,titleStr,...
    %             strcat(figurePrefixName,'_',tagDataType,'_',tagStatType),...
    %             figureSubfolderName,[])
    % 
    %     end
    % end % iDatatype
            
end % plotBcpMetrics

% *************************************************************************

function plotFlux2Dslices(figureSubfolderName,figurePrefixName,targetOutput,config)

    isAvgOrErr = 1; % 1 = average, 2 = error
    tagStatTypes = {'avg', 'err'};
    tagStatType = tagStatTypes{min(max(isAvgOrErr, 1), numel(tagStatTypes))};

    fluxMonthlySlicesMetadata = {...
        'POCflux_ex',   permute(squeeze(targetOutput.flux.monthly.dh(1,:,:,isAvgOrErr,1)),    [2,1]), 'POC flux at export depth',  'mg C m^{-2} d^{-1}',      0, 250;
        'POCflux_seq',  permute(squeeze(targetOutput.flux.monthly.dh(2,:,:,isAvgOrErr,1)),    [2,1]), 'POC flux at 1000 m',        'mg C m^{-2} d^{-1}',      0, 150;
        'PICflux_ex',   permute(squeeze(targetOutput.flux.monthly.dh(1,:,:,isAvgOrErr,2)),    [2,1]), 'PIC flux at export depth',  'mg CaCO3 m^{-2} d^{-1}',  0, 250;
        'PICflux_seq',  permute(squeeze(targetOutput.flux.monthly.dh(2,:,:,isAvgOrErr,2)),    [2,1]), 'PIC flux at 1000 m',        'mg CaCO3 m^{-2} d^{-1}',  0, 150;
        'bSiflux_ex',   permute(squeeze(targetOutput.flux.monthly.dh(1,:,:,isAvgOrErr,3)),    [2,1]), 'Opal flux at export depth', 'mg bSi m^{-2} d^{-1}',    0, 50;
        'bSiflux_seq',  permute(squeeze(targetOutput.flux.monthly.dh(2,:,:,isAvgOrErr,3)),    [2,1]), 'Opal flux at 1000 m',       'mg bSi m^{-2} d^{-1}',    0, 10;
        'POCtoPIC_ex',  permute(squeeze(targetOutput.fluxrat.monthly.dh(1,:,:,isAvgOrErr,1)), [2,1]), 'POC/PIC at export depth',   'mol C/mol C',             0, 30;
        'POCtoPIC_seq', permute(squeeze(targetOutput.fluxrat.monthly.dh(2,:,:,isAvgOrErr,1)), [2,1]), 'POC/PIC at 1000 m',         'mol C/mol C',             0, 30;
        'bSitoPIC_ex',  permute(squeeze(targetOutput.fluxrat.monthly.dh(1,:,:,isAvgOrErr,2)), [2,1]), 'bSi/PIC at export depth',   'mol Si/mol C',            0, 5;
        'bSitoPIC_seq', permute(squeeze(targetOutput.fluxrat.monthly.dh(2,:,:,isAvgOrErr,2)), [2,1]), 'bSi/PIC at 1000 m',         'mol Si/mol C',            0, 5;
    };

    fluxAnnualSlicesMetadata = {...
        'POCflux_ex',  permute(squeeze(targetOutput.flux.annual.dh(1,:,isAvgOrErr,1)), [2,1]), 'POC flux at export depth (mg C m^{-2} d^{-1})',     '', 0, 250;
        'POCflux_seq', permute(squeeze(targetOutput.flux.annual.dh(2,:,isAvgOrErr,1)), [2,1]), 'POC flux at 1000 m (mg C m^{-2} d^{-1})',           '', 0, 150;
        'PICflux_ex',  permute(squeeze(targetOutput.flux.annual.dh(1,:,isAvgOrErr,2)), [2,1]), 'PIC flux at export depth (mg CaCO3 m^{-2} d^{-1})', '', 0, 250;
        'PICflux_seq', permute(squeeze(targetOutput.flux.annual.dh(2,:,isAvgOrErr,2)), [2,1]), 'PIC flux at 1000 m (mg CaCO3 m^{-2} d^{-1})',       '', 0, 150;
        'bSiflux_ex',  permute(squeeze(targetOutput.flux.annual.dh(1,:,isAvgOrErr,3)), [2,1]), 'Opal flux at export depth (mg bSi m^{-2} d^{-1})',  '', 0, 50;
        'bSiflux_seq', permute(squeeze(targetOutput.flux.annual.dh(2,:,isAvgOrErr,3)), [2,1]), 'Opal flux at 1000 m (mg bSi m^{-2} d^{-1})',        '', 0, 10;
    };

    fluxRatiosAnnualSlicesMetadata = {...
        'POCtoPIC_ex',  permute(squeeze(targetOutput.fluxrat.annual.dh(1,:,isAvgOrErr,1)), [2,1]), 'POC/PIC at export depth (mol C/mol C)',  '', 0, 10;
        'POCtoPIC_seq', permute(squeeze(targetOutput.fluxrat.annual.dh(2,:,isAvgOrErr,1)), [2,1]), 'POC/PIC at 1000 m (mol C/mol C)',        '', 0, 10;
        'bSitoPIC_ex',  permute(squeeze(targetOutput.fluxrat.annual.dh(1,:,isAvgOrErr,2)), [2,1]), 'bSi/PIC at export depth (mol Si/mol C)', '', 0, 3;
        'bSitoPIC_seq', permute(squeeze(targetOutput.fluxrat.annual.dh(2,:,isAvgOrErr,2)), [2,1]), 'bSi/PIC at 1000 m (mol Si/mol C)',       '', 0, 3;
    };

    % Arrange into a structure the metadata
    dataTypes = struct( ...
        'metadata', {fluxMonthlySlicesMetadata, fluxAnnualSlicesMetadata, fluxRatiosAnnualSlicesMetadata}, ...
        'tag', {'output_global_monthly_flux', 'output_global_annual_flux', 'output_global_annual_flux_ratios'} ...
    );

    for iDatatype = 1:numel(dataTypes)
        datasetsMetadata = dataTypes(iDatatype).metadata;
        datasetsMetadata = unpackData(config,datasetsMetadata); % from location sequence structure to lat x lon array

        tagDataType = dataTypes(iDatatype).tag;

        % Initialise colour bar metadata
        cbString = cell(height(datasetsMetadata),1);
        caxisMin = zeros(height(datasetsMetadata),1);
        caxisMax = zeros(height(datasetsMetadata),1);
        for iDataset = 1:height(datasetsMetadata)
            caxisMin(iDataset) = datasetsMetadata{iDataset,5};
            caxisMax(iDataset) = datasetsMetadata{iDataset,6};
            cbString{iDataset} = datasetsMetadata{iDataset,4};
        end

        if iDatatype == 1 % monthly data (displayed as subplots in one figure)

            for iDataset = 1:height(datasetsMetadata)
                datasetName = lower(datasetsMetadata{iDataset,1});
                sgTitle = datasetsMetadata{iDataset,3};
                dataArray = datasetsMetadata{iDataset,2};
                plotOceanVariableMaps(dataArray,config.mappingProps.mapLons,...
                    config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
                    cbString{iDataset},caxisMin(iDataset),caxisMax(iDataset),...
                    [],config.mappingProps.labelMonths,...
                    strcat(figureSubfolderName,'_',figurePrefixName,'_',tagDataType,'_',tagStatType,'_',datasetName),...
                    figureSubfolderName,sgTitle)
            end

        else % annual data (group into subplots)

            groupDatasets = [];
            titleStr = cell(height(datasetsMetadata),1);
            for iDataset = 1:height(datasetsMetadata)
                dataArray = datasetsMetadata{iDataset,2};
                groupDatasets = cat(3,groupDatasets,dataArray); % append
                titleStr{iDataset} = datasetsMetadata{iDataset,3};
            end

            isCommonColourBar = false;
            plotOceanVariableMaps(groupDatasets,config.mappingProps.mapLons,...
                config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
                [],caxisMin,caxisMax,isCommonColourBar,titleStr,...
                strcat(figureSubfolderName,'_',figurePrefixName,'_',tagDataType,'_',tagStatType),...
                figureSubfolderName,[])

        end
    end % iDatatype

end % plotFlux2Dslices

% *************************************************************************

function plotAverageParticleProperties(figureSubfolderName,figurePrefixName,targetOutput,config)

    % Got rid of attribute #9, depth
    particleAttributesAnnualMetadata = {... % 3 sizes x 4 depth x nAtts x nLocs
        'pnum',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,1,:)), [3,1,2]),  'Particle number';
        'rho',   permute(squeeze(targetOutput.patt.annual.dh(2:3,:,2,:)), [3,1,2]),  'Density (g cm^{-3})';
        'velo',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,3,:)), [3,1,2]),  'Sinking speed (m d^{-1})';
        'stick', permute(squeeze(targetOutput.patt.annual.dh(2:3,:,4,:)), [3,1,2]),  'Stickiness';
        'por',   permute(squeeze(targetOutput.patt.annual.dh(2:3,:,5,:)), [3,1,2]),  'Porosity';
        'vol',   permute(squeeze(targetOutput.patt.annual.dh(2:3,:,6,:)), [3,1,2]),  'Volume (\mu^{3})';
        'frac',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,7,:)), [3,1,2]),  'Fractal dimension';
        'radpp', permute(squeeze(targetOutput.patt.annual.dh(2:3,:,8,:)), [3,1,2]),  'Radius primary particle (\mum)';
        'orgc',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,10,:)), [3,1,2]), 'POC content (pmol)';
        'tepc',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,11,:)), [3,1,2]), 'TEP-C content (pmol)';
        'opal',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,12,:)), [3,1,2]), 'Opal content (pmol)';
        'calc',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,13,:)), [3,1,2]), 'CaCO3 content (pmol)';
        'clay',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,14,:)), [3,1,2]), 'Clay content (pmol)';
        'size',  permute(squeeze(targetOutput.patt.annual.dh(2:3,:,15,:)), [3,1,2]), 'ESD (\mum)';
    };

    % From location sequence structure to lat x lon array
    particleAttributesAnnualMetadata = unpackData(config,particleAttributesAnnualMetadata); 
    
    % Define metadata for size categories and depth labels
    sizeCategories = ["Small particles (< 150 \mum ESD)", "Large particles (> 150 \mum ESD)"];
    depthLabels = ["Export depth", "500 m", "1000 m", "2000 m"];
    nSizeCategories = numel(sizeCategories);
    nDepthLevels = numel(depthLabels);
    nSubplots = nSizeCategories * nDepthLevels;

    for iAtt = 1:height(particleAttributesAnnualMetadata)
        attributeArray = particleAttributesAnnualMetadata{iAtt,2};
        attributeName = particleAttributesAnnualMetadata{iAtt,1};
        sgTitle = particleAttributesAnnualMetadata{iAtt,3};
        
        caxisMin = [];
        caxisMax = [];
        titleStr = cell(nSubplots,1);
        annotationStr = cell(nSubplots,1);
        groupDatasets = [];
        
        % Loop through size and depth indices
        for iDepth = 1:nDepthLevels
            for iSize = 1:nSizeCategories
                iSubplot = (iDepth - 1) * nSizeCategories + iSize;
                data = attributeArray(:,:,iSize,iDepth);
                titleStr{iSubplot} = sizeCategories(iSize);
                annotationStr{iSubplot} = depthLabels(iDepth);
                groupDatasets = cat(3,groupDatasets,squeeze(data)); % append
            end
        end 

%         if ((iAtt == 1 || iAtt == 7 || iAtt == 8 || iAtt == 9) && (iSubplot == 1 || iSubplot == 3 || iSubplot == 5 || iSubplot == 7))
%             mydata = log10(mydata);
%         end
        
        isCommonColourBar = false;
        plotOceanVariableMaps(groupDatasets,config.mappingProps.mapLons,...
            config.mappingProps.mapLats,config.mappingProps.myColourMapParticles,...
            annotationStr,caxisMin,caxisMax,isCommonColourBar,titleStr,...
            strcat(figureSubfolderName,'_',figurePrefixName,'_','output_global_annual_patt_',attributeName),...
            figureSubfolderName,sgTitle)

    end % iAtt
                
% geo_avgAttSc(geo_avgAttSc==0) = NaN;
% 
% for iAtt = 1:length(idxsAtts)
%         
%         if (iAtt == 1 && (iSubplot == 1 || iSubplot == 3 || iSubplot == 5 || iSubplot == 7))
%             set(gca,'ColorScale','log');
%         end
%         
%         cb = colorbar('Location','eastoutside');
%         if (iAtt == 1 && (iSubplot == 1 || iSubplot == 3 || iSubplot == 5 || iSubplot == 7))
%             cb.YTick = log10([1.5, 2, 10, 100]);  
%             cb.YTickLabel = {'1.5','2','10','100'};
%         elseif ((iAtt == 7 || iAtt == 8) && (iSubplot == 1 || iSubplot == 3 || iSubplot == 5 || iSubplot == 7))
%             cb.YTick = log10([1e-16, 1e-14, 1e-12, 1e-10]);  
%             cb.YTickLabel = {'1e-16','1e-14','1e-12','1e-10'};
%         elseif (iAtt == 9 && (iSubplot == 1 || iSubplot == 3 || iSubplot == 5 || iSubplot == 7))
%             cb.YTick = log10([1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6]);  
%             cb.YTickLabel = {'1e-16','1e-14','1e-12','1e-10','1e-8','1e-6'};
%         end   
% 
% end % iFigure

end % plotAverageParticleProperties

% *************************************************************************

function plotSmsAndAuxProperties(figureSubfolderName,figurePrefixName,targetOutput,config)
%%
    targetOutput.sms.annual.dh(targetOutput.sms.annual.dh==0) = NaN;
    targetOutput.aux.annual.dh(targetOutput.aux.annual.dh==0) = NaN;
    
    smsAnnualMetadata = {... % 33 idxs x 4 depths x nLocs
        'primprodorgc',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxPrimProdOrgC,:,:)), [2,1]),      'C taken up by phytoplankton';
        'primprodcaco3',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxPrimProdCaCO3,:,:)), [2,1]),     'Carbonate ion taken up by coccolithophores';
        'primprodopal',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxPrimProdOpal,:,:)), [2,1]),      'Silicic acid taken up by diatoms';
        'depoclay',         permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxDepoClay,:,:)), [2,1]),          'Aeolian dust deposition';
        'prodtepphyto',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxProdTepPhyto,:,:)), [2,1]),      'TEC released by phytoplankton during photosynthesis';
        'prodtepmicrob',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxProdTepMicrob,:,:)), [2,1]),     'TEC released by microbes';
        'phototepc',        permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxPhotoTepC,:,:)), [2,1]),         'DOC released due to photolysis';
        'zoosoluborgc',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooSolubOrgC,:,:)), [2,1]),      'DOC released due to zooplankton messy feeding';
        'zoosolubtepc',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooSolubTepC,:,:)), [2,1]),      'Dissolved TEC released due to zooplankton messy feeding';
        'zoosolubcaco3',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooSolubCaCO3,:,:)), [2,1]),     'CO3= released due to zooplankton messy feeding';
        'zoosolubopal',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooSolubOpal,:,:)), [2,1]),      'SiOH4 released due to zooplankton messy feeding';
        'zoosolubclay',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooSolubClay,:,:)), [2,1]),      'Dissolved clay released due to zooplankton messy feeding';
        'zooingestorgc',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooIngestOrgC,:,:)), [2,1]),     'Mesozoo. ingestion of POC';
        'zooingesttepc',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooIngestTepC,:,:)), [2,1]),     'Mesozoo. ingestion of TEC';
        'zooresporgc',      permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooRespOrgC,:,:)), [2,1]),       'C (CO2) released due to zooplankton respiration of POC';
        'zooresptepc',      permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooRespTepC,:,:)), [2,1]),       'C (CO2) released due to zooplankton respiration of TEC';
        'zooegestorgc',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooEgestOrgC,:,:)), [2,1]),      'Mesozoo. egestion of POC';
        'zooegesttepc',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooEgestTepC,:,:)), [2,1]),      'Mesozoo. egestion of TEC';
        'zoodissolcaco3',   permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooDissolCaCO3,:,:)), [2,1]),    'CO3= release due to biotic dissolution of CaCO3';
        'zooexcretorgc',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooExcretOrgC,:,:)), [2,1]),     'DOC release due to zooplankton excretion';
        'zooexcrettepc',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooExcretTepC,:,:)), [2,1]),     'Dissolved TEC release due to zooplankton excretion';
        'zoodeathorgc',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooDeathOrgC,:,:)), [2,1]),      'POC released as zooplankton dead bodies';
        'zoodeathcaco3',    permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooDeathCaCO3,:,:)), [2,1]),     'CaCO3 released as zooplankton dead bodies';
        'zoodeathopal',     permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxZooDeathOpal,:,:)), [2,1]),      'Opal released as zooplankton dead bodies';
        'microbresporgc',   permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxMicrobRespOrgC,:,:)), [2,1]),    'C (CO2) released due to microbial respiration of POC';
        'microbresptepc',   permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxMicrobRespTepC,:,:)), [2,1]),    'C (CO2) released due to microbial respiration of TEC';
        'microbsoluborgc',  permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxMicrobSolubOrgC,:,:)), [2,1]),   'DOC released due to microbial solubilisation';
        'microbsolubtepc',  permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxMicrobSolubTepC,:,:)), [2,1]),   'Dissolved TEC released due to microbial solubilisation';
        'microbsolubcaco3', permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxMicrobSolubCaCO3,:,:)), [2,1]),  'CO3= released due to abiotic solubilisation';
        'microbsolubopal',  permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxMicrobSolubOpal,:,:)), [2,1]),   'SiOH4 released due to abiotic solubilisation';
        'microbsolubclay',  permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxMicrobSolubClay,:,:)), [2,1]),   'Dissolved clay released due to abiotic solubilisation';
        'dissolcaco3',      permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxDissolCaCO3,:,:)), [2,1]),       'CO3= released due to abiotic dissolution of CaCO3';
        'dissolopal',       permute(squeeze(targetOutput.sms.annual.dh(config.sms.idxDissolOpal,:,:)), [2,1]),        'SiOH4 released due to abiotic dissolution of opal';
    };

    auxAnnualMetadata = {... % 32 idxs x 4 depths x nLocs
        'zoobiomass',                   permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxZooBiomass,:,:)), [2,1]),                'Mesozoo. biomass',                                   'mol C m^{-3}';
        'zoonumber',                    permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxZooNumber,:,:)), [2,1]),                 'No. mesozoo. individuals',                           '# ind. m^{-3}';
        'zoodeadbiomass',               permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxZooDeadBiomass,:,:)), [2,1]),            'Mesozoo. dead biomass',                              'mol C m^{-3}';
        'diatbiomass',                  squeeze(targetOutput.aux.annual.dh(config.aux.idxDiatBiomass,1,:)),                               'Diatom biomass',                                     'mol C m^{-3}';
        'flagelbiomass',                squeeze(targetOutput.aux.annual.dh(config.aux.idxFlagelBiomass,1,:)),                             'Large non-mineralising phytoplankton biomass',       'mol C m^{-3}';
        'coccobiomass',                 squeeze(targetOutput.aux.annual.dh(config.aux.idxCoccoBiomass,1,:)),                              'Coccolithophore biomass',                            'mol C m^{-3}';
        'picobiomass',                  squeeze(targetOutput.aux.annual.dh(config.aux.idxPicoBiomass,1,:)),                               'Picophytoplankton',                                  'mol C m^{-3}';
        'freshdiatcellquota',           squeeze(targetOutput.aux.annual.dh(config.aux.idxFreshDiatCellQuota,1,:)),                        'Fresh diatoms cell quota',                           'mol C cell^{-1}';
        'freshflagelcellquota',         squeeze(targetOutput.aux.annual.dh(config.aux.idxFreshFlagelCellQuota,1,:)),                      'Fresh large non-mineralising phytoplankton biomass', 'mol C cell^{-1}';
        'freshcoccocellquota',          squeeze(targetOutput.aux.annual.dh(config.aux.idxFreshCoccoCellQuota,1,:)),                       'Fresh coccolithophore cell quota',                   'mol C cell^{-1}';
        'freshpicocellquota',           squeeze(targetOutput.aux.annual.dh(config.aux.idxFreshPicoCellQuota,1,:)),                        'Fresh picophytoplankton cell quota',                 'mol C cell^{-1}';
        'zoospecresprate',              permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxZooSpecRespRate,:,:)), [2,1]),           'Mesozoo. specific respiration rate',                 's^{-1}';
        'microbspecresprate',           permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxMicrobSpecRespRate,:,:)), [2,1]),        'Microbial specific respiration rate',                's^{-1}';
        'nightdvmupperbound',           squeeze(targetOutput.aux.annual.dh(config.aux.idxNightDvmUpperBound,1,:)),                        'Night DVM upper boundary',                           'm';
        'nightdvmlowerbound',           squeeze(targetOutput.aux.annual.dh(config.aux.idxNightDvmLowerBound,1,:)),                        'Night DVM lower boundary',                           'm';
        'daydvmupperbound',             squeeze(targetOutput.aux.annual.dh(config.aux.idxDayDvmUpperBound,1,:)),                          'Day DVM upper boundary',                             'm';
        'daydvmlowerbound',             squeeze(targetOutput.aux.annual.dh(config.aux.idxDayDvmLowerBound,1,:)),                          'Day DVM lower boundary',                             'm';
        'euphoticdepth',                squeeze(targetOutput.aux.annual.dh(config.aux.idxEuphoticDepth,1,:)),                             'Euphotic depth',                                     'm';
        'collisionkernel',              permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxCollisionKernel,:,:)), [2,1]),           'Collision kernel',                                   'm^{-3} s^{-1}';
        'browniankernel',               permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxBrownianKernel,:,:)), [2,1]),            'Brownian kernel',                                    'm^{-3} s^{-1}';
        'shearkernel',                  permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxShearKernel,:,:)), [2,1]),               'Shear rate kernel',                                  'm^{-3} s^{-1}';
        'settlingkernel',               permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxSettlingKernel,:,:)), [2,1]),            'Settling kernel',                                    'm^{-3} s^{-1}';
        'numclusterpairsevalcoagu',     permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumClusterPairsEvalCoagu,:,:)), [2,1]),  'No. cluster pairs evaluated for coagulation',        '';
        'numpxcmore',                   permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumPxCmore,:,:)), [2,1]),                'No. particles per cluster in cluster with more part.','';
        'numpxcless',                   permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumPxCless,:,:)), [2,1]),                'No. particles per cluster in cluster with less part.','';
        'coagulationsuccess',           permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxCoagulationSuccess,:,:)), [2,1]),        'Coagulation success',                                '';
        'coagulationprobability',       permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxCoagulationProbability,:,:)), [2,1]),    'Coagulation probability',                            '';
        'numparticlesevalencounter',    permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumParticlesEvalEncounter,:,:)), [2,1]), 'No. particles evaluated for encounter',              '';
        'numparticlesencountered',      permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumParticlesEncountered,:,:)), [2,1]),   'No. particles encountered by mesozoo.',              '';
        'numparticlesfragmentedzoo',    permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumParticlesFragmentedZoo,:,:)), [2,1]), 'No. particles fragmented by mesozoo.',               '';
        'numparticlesingestedzoo',      permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumParticlesIngestedZoo,:,:)), [2,1]),   'No. particles ingested by mesozoo.',                 '';
        'numparticlesomittedzoo',       permute(squeeze(targetOutput.aux.annual.dh(config.aux.idxNumParticlesOmittedZoo,:,:)), [2,1]),    'No. particles omitted by mesozoo.',                  '';
    };

    % From location sequence structure to lat x lon array
    smsAnnualMetadata = unpackData(config,smsAnnualMetadata);
    auxAnnualMetadata = unpackData(config,auxAnnualMetadata);
    
    % Define metadata for depth labels
    depthLabels = ["Export depth", "500 m", "1000 m", "2000 m"];
    nDepthLevels = numel(depthLabels);

    % Plot SMS terms
    for iSms = 1:height(smsAnnualMetadata)
        smsArray = smsAnnualMetadata{iSms,2};
        smsName = smsAnnualMetadata{iSms,1};
        sgTitle = smsAnnualMetadata{iSms,3};
        cbString = 'mol m^{-3} d^{-1}';
        caxisMin = min(smsArray,[],'all'); 
        caxisMax = max(smsArray,[],'all'); 
        
        if isnan(caxisMin) % means that SMS term is not being used any more by current SLAMS version
            continue;
        end

        titleStr = cell(nDepthLevels,1);
        groupDatasets = [];
        
        % Loop through depth indices
        for iDepth = 1:nDepthLevels
            data = smsArray(:,:,iDepth);
            titleStr{iDepth} = depthLabels(iDepth);
            groupDatasets = cat(3,groupDatasets,squeeze(data)); % append
        end 

        isCommonColourBar = true;
        plotOceanVariableMaps(groupDatasets,config.mappingProps.mapLons,...
            config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
            cbString,caxisMin,caxisMax,isCommonColourBar,titleStr,...
            strcat('output_global_annual_sms_',smsName),sgTitle)

    end % iSms

    % Plot auxialliary terms 
    singleLayerTerms = [config.aux.idxDiatBiomass, config.aux.idxFlagelBiomass, ...
                        config.aux.idxCoccoBiomass, config.aux.idxPicoBiomass, ...
                        config.aux.idxFreshDiatCellQuota, config.aux.idxFreshFlagelCellQuota, ...
                        config.aux.idxFreshCoccoCellQuota, config.aux.idxFreshPicoCellQuota, ...
                        config.aux.idxNightDvmUpperBound, config.aux.idxNightDvmLowerBound, ...
                        config.aux.idxDayDvmUpperBound, config.aux.idxDayDvmLowerBound, ...
                        config.aux.idxEuphoticDepth];
    depthTerms = [config.aux.idxNightDvmUpperBound, config.aux.idxNightDvmLowerBound, ...
                  config.aux.idxDayDvmUpperBound, config.aux.idxDayDvmLowerBound, ...
                  config.aux.idxEuphoticDepth];
                
    for iAux = 1:height(auxAnnualMetadata)
        auxArray = auxAnnualMetadata{iAux,2};
        auxName = auxAnnualMetadata{iAux,1};
        cbString = auxAnnualMetadata{iAux,4};
        caxisMin = min(auxArray,[],'all'); 
        caxisMax = max(auxArray,[],'all'); 
        
        if isnan(caxisMin) % means that aux. term is not being used any more by current SLAMS version
            continue;
        end
        
        if ismember(iAux,singleLayerTerms)
            
            if ismember(iAux,depthTerms)
                % Round up to avoid scaling issues
                auxArray(auxArray < 0.1) = 0.1;
            end
            
            titleStr = auxAnnualMetadata{iAux,3};
            plotOceanVariableMaps(auxArray,config.mappingProps.mapLons,...
                config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
                cbString,caxisMin,caxisMax,[],titleStr,...
                strcat('output_global_annual_aux_',auxName),[])
            
        else
            
            sgTitle = auxAnnualMetadata{iAux,3};        
            titleStr = cell(nDepthLevels,1);
            groupDatasets = [];

            % Loop through depth indices
            for iDepth = 1:nDepthLevels
                data = auxArray(:,:,iDepth);
                titleStr{iDepth} = depthLabels(iDepth);
                groupDatasets = cat(3,groupDatasets,squeeze(data)); % append
            end 

            isCommonColourBar = true;
            plotOceanVariableMaps(groupDatasets,config.mappingProps.mapLons,...
                config.mappingProps.mapLats,config.mappingProps.myColourMapOceanVars,...
                cbString,caxisMin,caxisMax,isCommonColourBar,titleStr,...
                strcat('output_global_annual_aux_',auxName),sgTitle)
            
        end

    end % iAux
%%    
end % plotSmsAndAuxProperties

% *************************************************************************

function plotContourPlotsForcingVsMetrics(figureSubfolderName,figurePrefixName,model,forcingDataGeo,targetOutput,config)

    % Surface properties (good as they are)
    geoAnnualNpp = forcingDataGeo.annualNpp.*(3600*24*config.molarMassCarbon*1e3); % mol C m-2 s-1 --> mg C m-2 d-1 
    annualNpp = geoAnnualNpp(:); % from 2D array to single column vector
    annualChla = forcingDataGeo.annualChla(:); % from 2D array to single column vector
    annualPar0 = forcingDataGeo.annualPar0(:); % from 2D array to single column vector

    % Zooplankton: integrate over depth using trapezoidal integration
    geoAnnualMesozooDepthIntegrated = integrateOverDepth(forcingDataGeo.annualMesozoo,config); % mg C m-3 --> mg C m-2
    annualMesozooDepthIntegrated = geoAnnualMesozooDepthIntegrated(:); % from 2D array to single column vector

    % Average surface water-column nitrate and temperature
    geoAnnualMeanSurfTemp = NaN(config.nLats,config.nLons);
    geoAnnualMeanSurfNit = NaN(config.nLats,config.nLons);
    for iLat = 1:config.nLats
        for iLon = 1:config.nLons
            [~,iDepth] = min(abs(config.geoDepths(iLat,iLon,:) - 200)); % like for zooplankton, we'll use depth limit of 200 m
            if ~isnan(iDepth)
                geoAnnualMeanSurfTemp(iLat,iLon) = mean(forcingDataGeo.annualTemp(iLat,iLon,1:iDepth),'omitnan');
                geoAnnualMeanSurfNit(iLat,iLon) = mean(forcingDataGeo.annualNit(iLat,iLon,1:iDepth),'omitnan');
            end
        end
    end

    annualMeanSurfTemp = geoAnnualMeanSurfTemp(:); % from 2D array to single column vector
    annualMeanSurfNit = geoAnnualMeanSurfNit(:); % from 2D array to single column vector

    % Put them together
    forcingArrays = {annualNpp, annualChla, annualPar0, annualMeanSurfTemp, annualMeanSurfNit, annualMesozooDepthIntegrated};

    % Relative phytoplankton biomass
    biomassIndices = [config.aux.idxDiatBiomass, config.aux.idxFlagelBiomass, config.aux.idxCoccoBiomass, config.aux.idxPicoBiomass];
    phytoBiomass = squeeze(targetOutput.aux.annual.dh(biomassIndices,1,:))';  % 4 x nLocs
    totalPhytoBiomass = sum(phytoBiomass,2); % sum across PFTs
    pftArray = phytoBiomass ./ totalPhytoBiomass; % element-wise division

    % % Option B
    % pftArray = geo.annualProbPft;
    % pftArray = relatAbundances(:); % one column vector

    metricsMetadata = {...
        'peeff',    model.peeff.annual,     'Particle export efficiency',       [0, 0.3];
        'teff',     model.teff.annual,      'Transfer efficiency',              [0, 0.5];
        'martinb',  model.martinb.annual,   'Martin b',                         [-4, 4];
        'zstar',    model.zstar.annual,     'Remineralisation length scale',    [0, 500];
        'xi',       model.xi.annual,        'Slope particle size distribution', [-7, 7];
    };

    nSubplots = length(forcingArrays)*size(pftArray,2); % should be 6 x 4
    contourLevels = 0:0.05:0.8; % Specify the levels you want

    for iFigure = 1:5 % metrics to plot

        zvar = metricsMetadata{iFigure,2};
        cbString = metricsMetadata{iFigure,3};
        figureName = metricsMetadata{iFigure,1};
        caxisMin = metricsMetadata{iFigure,4}(1);
        caxisMax = metricsMetadata{iFigure,4}(2);

        figure()
        set(gcf,'Units','Normalized','Position',[0.0 0.0 0.55 0.80],'Color','w') 
        haxis = zeros(nSubplots,1);

        for iSubplot = 1:nSubplots

            haxis(iSubplot) = subaxis(6,4,iSubplot,'Spacing',0.011,'Padding',0.011,'Margin',0.07);
            ax(iSubplot).pos = get(haxis(iSubplot),'Position');

            ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.02; 
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

            if (iSubplot >= 1 && iSubplot <= 4)
                yvar = forcingArrays{1};
                ylabelStr = 'NPP (mg C m^{-2} d^{-1})';
            elseif (iSubplot >=5 && iSubplot <= 8)
                yvar = forcingArrays{2};
                ylabelStr = 'Chla (mg m^{-3})';
            elseif (iSubplot >=9 && iSubplot <= 12)
                yvar = forcingArrays{3};
                ylabelStr = 'PAR0 (W m^{-2})';
            elseif (iSubplot >=13 && iSubplot <= 16)
                yvar = forcingArrays{4};
                ylabelStr = 'Temperature (ºC)';
            elseif (iSubplot >=17 && iSubplot <= 20)
                yvar = forcingArrays{5};
                ylabelStr = 'Nitrate (mmol m^{-3})';
            elseif (iSubplot >=21 && iSubplot <= 24)
                yvar = forcingArrays{6};
                ylabelStr = 'Mesozoo. (mg C m^{-2})';
            end

            if (iSubplot == 1 || iSubplot == 5 || iSubplot == 9 || iSubplot == 13 || iSubplot == 17 || iSubplot == 21)
                xvar = 100.*pftArray(:,1); 
                xlabelStr = 'Diatom rel. abund. (%)';
                xMin = 16;
                xMax = 55;
            elseif (iSubplot == 2 || iSubplot == 6 || iSubplot == 10 || iSubplot == 14 || iSubplot == 18 || iSubplot == 22)
                xvar = 100.*pftArray(:,2);
                xlabelStr = 'Flagellate rel. abund. (%)';
                xMin = 7;
                xMax = 20.4;
            elseif (iSubplot == 3 || iSubplot == 7 || iSubplot == 11 || iSubplot == 15 || iSubplot == 19 || iSubplot == 23)
                xvar = 100.*pftArray(:,3); 
                xlabelStr = 'Coccolithophore rel. abund. (%)';
                xMin = 9;
                xMax = 33;
            elseif (iSubplot == 4 || iSubplot == 8 || iSubplot == 12 || iSubplot == 16 || iSubplot == 20 || iSubplot == 24)
                xvar = 100.*pftArray(:,4); 
                xlabelStr = 'Picophytoplankton rel. abund. (%)';
                xMin = 16;
                xMax = 54;
            end

            % Sort data and remove NaNs
            [xvar_sort, sortIdx] = sort(xvar);
            yvar_sort            = yvar(sortIdx);
            zvar_sort            = zvar(sortIdx);

        %     idxFirstZeroXvar   = find(xvar_sort == 0, 1, 'last');
        %     xvar_sort          = xvar_sort(idxFirstZeroXvar+1:end);
        %     yvar_sort          = yvar_sort(idxFirstZeroXvar+1:end);
        %     zvar_sort          = zvar_sort(idxFirstZeroXvar+1:end);

            idxFirstNanXvar = find(isnan(xvar_sort), 1, 'first');
            if ~isempty(idxFirstNanXvar)
                xvar_sort = xvar_sort(1:idxFirstNanXvar-1);
                yvar_sort = yvar_sort(1:idxFirstNanXvar-1);
                zvar_sort = zvar_sort(1:idxFirstNanXvar-1);
                yvar_sort(isnan(yvar_sort)) = 0;
                zvar_sort(isnan(zvar_sort)) = 0;
            end

            xvar_sort(isnan(xvar_sort)) = 0;
            yvar_sort(isnan(yvar_sort)) = 0;
            zvar_sort(isnan(zvar_sort)) = 0;

            minxvar = min(xvar_sort(:));
            maxxvar = max(xvar_sort(:));
            minyvar  = min(yvar_sort(:));
            maxyvar  = max(yvar_sort(:));

            % Interpolation grid
            xv = linspace(minxvar, maxxvar, 150);
            yv = linspace(minyvar, maxyvar, 150);
            [Xm,Ym] = ndgrid(xv, yv);
            Zm = griddata(xvar_sort, yvar_sort, zvar_sort, Xm, Ym);

            % Plot contour
            contourf(haxis(iSubplot),Xm,Ym,Zm,contourLevels,'LineStyle','none')
            colormap(jet(length(contourLevels)-1))
            caxis([caxisMin caxisMax])
            shading flat
            xlim([xMin xMax])

            if (iSubplot >= 1 && iSubplot <= 4)
                ylim([100 1000])
            elseif (iSubplot >=5 && iSubplot <= 8)
                ylim([0 3])
            end

            if (iSubplot >= 21 && iSubplot <= 24)
                xh = xlabel(xlabelStr);
                xh.Position(2) = xh.Position(2) - 0.30;
            else
                xlabel([]);
            end

            if (iSubplot == 1 || iSubplot == 5 || iSubplot == 9 || iSubplot == 13 || iSubplot == 17 || iSubplot == 21)
                yh = ylabel(ylabelStr);
                yh.Position(1) = yh.Position(1) - 1;
            else
                ylabel([]); 
            end
        %     th.Position(2) = th.Position(2) + 30;

            if (iSubplot == 4)
                cb = colorbar('Location','eastoutside');
                cb.Position(1) = cb.Position(1) + 0.07;
                cb.Position(2) = cb.Position(2) - 0.73;
                cb.Position(3) = 0.015; % WIDTH
                cb.Position(4) = 0.845; % LENGTH
                cb.Label.String = cbString;
                cb.FontSize = 12;
            end

            grid on

            set(haxis(iSubplot),'FontSize',10)

        end

        exportgraphics(gcf,fullfile('.','figures',...
            strcat('contour_pfts_',figureName,'.png')),'Resolution',600) 

    end % iFigure

end % plotContourPlotsForcingVsMetrics

% *************************************************************************

function data = unpackData(config,data)

    % Loop through each entry in the metadata
    for iField = 1:size(data,1)

        % Extract data field
        dataField = data{iField,2}; 

        % Loop through each location
        for iLoc = 1:size(dataField,1)

            % Extract location indices
            iLat = config.iyBb(iLoc);
            iLon = config.ixBb(iLoc);

            % Get dimensions of the data, excluding location
            sizeDataDims = size(dataField);
            sizeDimsExcludingLoc = sizeDataDims(2:end);

            % Build dynamic index structure for data extraction
            idxsInDim = arrayfun(@(dimSize) 1:dimSize, sizeDimsExcludingLoc, 'UniformOutput', false);

            % Prepend location index
            idxsInDim = [{iLoc}, idxsInDim];

            % Extract the data at this location
            dataLoc = dataField(idxsInDim{:});

            % Build dynamic index structure
            seqIndices = [{iLat},{iLon},repmat({':'},1,length(sizeDimsExcludingLoc))];
            
            if iLoc == 1
                newDataArray = NaN([config.nLats, config.nLons, sizeDimsExcludingLoc]);
            end

            % Assign the modified data **back** into the original array
            newDataArray(seqIndices{:}) = dataLoc;
            
        end % iLoc
        
        data{iField,2} = newDataArray;

    end % iField

end % unpackData

% *************************************************************************

function [qZeuMonthly,qZeuAnnual,qNppMonthly] = extractOutputSupportingData(...
    config,filenameInputZeu,filenameInputNpp,filenameInputMask,logID)

    % Load the global-ocean euphotic layer depth product calculated from
    % chla from OC-CCI and MLD from IFREMER
    load(fullfile('.','data','interim',filenameInputZeu),'zeu','zeu_lat','zeu_lon')
 
    % Load NPP from BICEP
    load(fullfile('data','raw',filenameInputNpp),'npp_avg','npp_err','npp_lat','npp_lon')
    
    % Load mask
    load(fullfile('data','raw',filenameInputMask),'mask','mask_lat','mask_lon')
    
    % Mask NPP
    nppAvgMasked = manageGapsInOceanDataIceMasked(npp_avg,mask,npp_lat,npp_lon,...
        mask_lat,mask_lon,(1:12),logID);
    nppErrMasked = manageGapsInOceanDataIceMasked(npp_err,mask,npp_lat,npp_lon,...
        mask_lat,mask_lon,(1:12),logID);
    clear mask
    
    % Original data grids for zeu and NPP
    [Xz,Yz,Tz] = ndgrid(zeu_lat,zeu_lon,(1:12)');
    [Xn,Yn,Tn] = ndgrid(npp_lat,npp_lon,(1:12)');

    % Interpolants
    Fz = griddedInterpolant(Xz, Yz, Tz, zeu, 'linear'); 
    Fna = griddedInterpolant(Xn, Yn, Tn, nppAvgMasked, 'linear'); 
    Fne = griddedInterpolant(Xn, Yn, Tn, nppErrMasked, 'linear');

    % Extract data for the study locations
    qZeuMonthly = NaN(12,config.nLocs,2); % 3rd dim: 1=value, 2=err
    qZeuAnnual  = NaN(config.nLocs,2); % 2nd dim: 1=value, 2=err
    qNppMonthly = NaN(12,config.nLocs,2); 

    for iLoc = 1:config.nLocs
        iLat = config.iyBb(iLoc);
        iLon = config.ixBb(iLoc);
        
        [qX,qY,qT] = ndgrid(config.lats(iLat),config.lons(iLon),(1:12)');
        
        % Value at location
        zeuLocal = squeeze(Fz(qX,qY,qT));
        nppAvgLocal = squeeze(Fna(qX,qY,qT));
        nppErrLocal = squeeze(Fne(qX,qY,qT));

        % Fill NaN gaps strategically for zeu (NPP is already good)
        zeuLocal = manageGapsInOceanData(zeuLocal,[],...
            (1:12),'false',Fz,config.lats(iLat),config.lons(iLon),...
            config.lonIncrement,config.latIncrement,config.lats(end),...
            config.maxIterForNeighbourhoodExpansion,logID);

        % Add uncertainty
        qZeuMonthly(:,iLoc,1) = zeuLocal; 
        qZeuMonthly(:,iLoc,2) = config.errorFractionZeu.*zeuLocal;
        qZeuAnnual(iLoc,1) = mean(zeuLocal,'omitnan');
        qZeuAnnual(iLoc,2) = mean(config.errorFractionZeu.*zeuLocal,'omitnan');
        
        qNppMonthly(:,iLoc,1) = nppAvgLocal;
        qNppMonthly(:,iLoc,2) = nppErrLocal;
    end

end % extractOutputSupportingData

% *************************************************************************

function [metricAnnual,metricMonthly] = calculateEfficiencyMetricsWithoutErrorPropagation(...
    arrayNumerator,arrayDenominator)

nLocs = size(arrayNumerator,2);
metricMonthly = NaN(nLocs,12); 
metricAnnual = NaN(nLocs,1);  

% Calculate the quotient for valid entries
valid = ~isnan(arrayNumerator) & ~isnan(arrayDenominator) & (arrayDenominator ~= 0);
ratio = NaN(size(arrayNumerator));
ratio(valid) = arrayNumerator(valid)./arrayDenominator(valid);

metricMonthly(:,:) = ratio';
metricAnnual(:)    = squeeze(nanmean(ratio,1)); % mean over months

end % calculateEfficiencyMetricsWithoutErrorPropagation

% *************************************************************************

