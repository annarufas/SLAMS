
% ======================================================================= %
%                                                                         %
% This script plots all relevant plots for parameter sensitivity analysis.%
%                                                                         %
%   Version 1.0 - Completed 29 May 2025                                   %
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

% Name of test cases
testrunDir = {'LOCALTS6';... % default run - 1
              'LOCALTS6_Si2C_diat_LOW'; 'LOCALTS6_Si2C_diat_UPP';... % 2-3
              'LOCALTS6_Calc2C_cocco_max_LOW'; 'LOCALTS6_Calc2C_cocco_max_UPP';... % 4-5
              'LOCALTS6_phyto_exudation_frac_LOW'; 'LOCALTS6_phyto_exudation_frac_UPP';... % 6-7
              'LOCALTS6_zoo_distrib_slope_LOW'; 'LOCALTS6_zoo_distrib_slope_UPP';... % 8-9
              'LOCALTS6_detection_radius_factor_mesozoo_LOW'; 'LOCALTS6_detection_radius_factor_mesozoo_UPP';... % 10-11
              'LOCALTS6_agg_to_zoo_size_ratio_LOW'; 'LOCALTS6_agg_to_zoo_size_ratio_UPP';... % 12-13
              'LOCALTS6_frac_OM_zoo_ingestion_surf_LOW'; 'LOCALTS6_frac_OM_zoo_ingestion_surf_UPP';... % 14-15
              'LOCALTS6_q10_microb_LOW'; 'LOCALTS6_q10_microb_UPP';... % 16-17
              'LOCALTS6_resp_rate_poc_max_0deg_microb_LOW'; 'LOCALTS6_resp_rate_poc_max_0deg_microb_UPP';... % 18-19
              'LOCALTS6_resp_rate_tepc_max_0deg_microb_LOW'; 'LOCALTS6_resp_rate_tepc_max_0deg_microb_UPP';... % 20-21
              'LOCALTS6_k_O2_resp_LOW'; 'LOCALTS6_k_O2_resp_UPP';... % 22-23
              'LOCALTS6_solub_rate_poc_LOW'; 'LOCALTS6_solub_rate_poc_UPP';... % 24-25
              'LOCALTS6_solub_rate_tepc_LOW'; 'LOCALTS6_solub_rate_tepc_UPP';... % 26-27
              'LOCALTS6_dissol_rate_calc_LOW'; 'LOCALTS6_dissol_rate_calc_UPP';... % 28-29
              'LOCALTS6_dissol_rate_opal_0deg_LOW'; 'LOCALTS6_dissol_rate_opal_0deg_UPP';... % 30-31
              'LOCALTS6_q10_bSi_LOW'; 'LOCALTS6_q10_bSi_UPP';... % 32-33
              'LOCALTS6_photodegradation_rate_tepc_LOW'; 'LOCALTS6_photodegradation_rate_tepc_UPP'}; % 34-35

paramMetadata = {...
    '$$\theta_{\mathrm{bSi:C,diat}}$$',                   0.13,   0.01,   0.50;   % Diatom cellular opal:POC
    '$$\theta_{\mathrm{calc:C,cocco}}^{\mathrm{max}}$$',  0.50,   0.05,   2.0;    % Max. coccolithophore cellular PIC:POC
    '$$f_{\mathrm{TEP}}$$',                               0.40,   0.02,   0.60;   % Fraction of a phytoplankton cell exuded as TEP matter
    '$$\xi_{\mathrm{z}}$$',                               1.0,    0.40,   1.0;    % Slope of the mesozooplankton normalised biomass size spectra
    '$$R_{\mathrm{z}}$$',                                 2.64,   1.0,    5.0;    % Mesozooplankton detection radius factor
    '$$\theta_{\mathrm{agg:zoo}}$$',                      1.0,    0.20,   4.0;    % Ratio aggregate to zooplankter
    '$$f_{\mathrm{OM}}^{\mathrm{ing,surf}}$$',            0.20,   0.10,   0.80;   % Min. fraction of OM in a surface-ocean particle to be deemed appetitive
    '$$Q_{\mathrm{10,m}}$$',                              2.2,    1.5,    3.0;    % Microbial Q10
    '$$\lambda_{\mathrm{m,C_{POC}}}^{\mathrm{max,0}}$$',  0.01,   0.001,  0.20;   % Max. microbial POC-specific respiration rate at 0ºC
    '$$\lambda_{\mathrm{m,C_{TEP}}}^{\mathrm{max,0}}$$',  0.01,   0.001,  0.10;   % Max. microbial TEP-C-specific respiration rate at 0ºC
    '$$k_{\mathrm{O2}}$$',                                0.50,   0.0005, 1.0;    % O2 half-saturation constant for respiration
    '$$\omega_{\mathrm{C_{POC}}}$$',                      0.005,  0.001,  0.10;   % Microbial solubilisation rate for Corg
    '$$\omega_{\mathrm{C_{TEP}}}$$',                      0.010,  0.001,  0.10;   % Microbial solubilisation rate for TEP-C
    '$$\kappa_{\mathrm{calc}}$$',                         2.0,    1.0,    7.0;    % Calcite dissolution rate
    '$$\kappa_{\mathrm{bSi}}^{\mathrm{max,0}}$$',         0.10,   0.0005, 0.50;   % Max. opal dissolution rate at 0ºC
    '$$Q_{\mathrm{10,bSi}}$$',                            2.3,    1.5,    3.5;    % Opal dissolution rate coefficient
    '$$\phi_{\mathrm{photo}}$$',                          0.30,   0.10,   0.70;   % Photodegradation rate of TEP
}; 

figureSubfolder = 'test_param_sensitivity';

% Choices to be applied
choiceTypeGridDomain = 2; % 1=global, 2=local, 3=global with zoom in factor

% Filename declarations
filenameSlamsOutput           = 'runsoutput.mat';
filenameSlamsRunGrid          = 'grid_run.mat';
filenameSlamsNumDepthLayers   = 'waterColNumDepthLayers.txt';

filenameTimeseriesInformation = 'timeseries_station_information_slams.mat';
filenameObsPnumUvp5           = 'pnum_16sc_compilation_slams.mat';
filenameObsPocFlux            = 'pocflux_compilation_slams.mat';
filenameObsPicAndBsiFlux      = 'picandbsiflux_compilation_slams.mat';
filenameBicepPft              = 'pft_bicep.mat';
filenameBicepExportFluxDunne  = 'ef_dunne_bicep.mat';
filenameBicepExportFluxHenson = 'ef_henson_bicep.mat';
filenameBicepExportFluxLi     = 'ef_li_bicep.mat';
filenameZeu                   = 'zeu_calculated_chlaoccci_mldifremer_pointonepercentpar0.mat';

% Path declarations
fullpathTestDir             = fullfile('.','tests',testrunDir);
fullpathRawDataDir          = fullfile('.','data','raw');
fullpathModelInputDataDir   = fullfile(fullpathTestDir,'modelinputdata');
fullpathModelRunsDir        = fullfile(fullpathTestDir,'modelruns');

% Rearrange model locations to match observational order
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'}; % as in config.gridLats, config.gridLons
desiredModLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'}; % as in config.gridLats, config.gridLons
[~,reorderModLocIdx] = ismember(desiredModLocationOrder,currentModLocationOrder); % get reordering indices
    
% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - LOAD SELECTED MODEL OUTPUT AND EXTRACT DATA AT RELEVANT
% DEPTHS
% -------------------------------------------------------------------------

for iTest = 1:length(testrunDir)
    currentDir = testrunDir{iTest};
    fullpathTestDir           = fullfile('.','tests',currentDir);
    fullpathModelRunsDir      = fullfile(fullpathTestDir,'modelruns');
    fullpathModelInputDataDir = fullfile(fullpathTestDir,'modelinputdata');
    
    % Load configuration parameters
    if iTest == 1
        config = loadModelConfigurationParameters(fullpathModelInputDataDir,...
            filenameSlamsRunGrid,filenameSlamsNumDepthLayers,choiceTypeGridDomain);
    end

    % Load modelled flux data and format it to pair it with observations
    [~,modFluxMonthlyDh,modFluxMonthlyDhDepths,~,~,~,obsDataMonthly,obsErrorMonthly,...
     ~,~,bicepExFluxMonthly,~,~,~] =...
        loadAndFormatModelledFluxDataAndObservations(choiceTypeGridDomain,config,...
            fullpathModelRunsDir,filenameSlamsOutput,filenameObsPocFlux,...
            filenameObsPicAndBsiFlux,filenameTimeseriesInformation,...
            filenameBicepExportFluxDunne,filenameBicepExportFluxHenson,...
            filenameBicepExportFluxLi,filenameZeu);

    % Load particle number data
    [modPnumByScMonthlyProfile,~,~,~,uvpPnumMonthlyProfile,~,uvpPnumMonthlyTargetValues,...
        uvpPnumMonthlyTargetDepths,~,~] =...
        loadAndFormatModelledParticleNumberDataAndObservations(choiceTypeGridDomain,config,...
            filenameTimeseriesInformation,filenameObsPnumUvp5,fullpathModelRunsDir,filenameSlamsOutput);

    % Sum particle number data across size classes
    modTotPnumMonthlyProfile = squeeze(sum(modPnumByScMonthlyProfile,1,'omitnan'));
    uvpTotPnumMonthlyTarget = squeeze(sum(uvpPnumMonthlyTargetValues,1,'omitnan')); % the depths in uvpPnumMonthlyTargetValues match those in modPnumMonthlyProfile 

    % Need to transform modTotPnumMonthlyProfile and uvpTotPnumMonthlyTarget 
    % into something more practical (i.e., with less depths). For that,
    % let's target closest depths to modFluxMonthlyDhDepths
    nTargetDepths = size(modFluxMonthlyDh,1);
    modTotPnumMonthlyDh = NaN(nTargetDepths,12,config.nLocs);
    uvpTotPnumMonthlyDh = NaN(nTargetDepths,12,config.nLocs); 

    for iDh = 1:nTargetDepths
        for iMonth = 1:12
            for iLoc = 1:config.nLocs

                % Find the closest depth to modelledTargetDepths
                uvpDepths = squeeze(uvpPnumMonthlyTargetDepths(1,:,iMonth,iLoc));
                modelledTargetDepths = squeeze(modFluxMonthlyDhDepths(iDh,iMonth,iLoc,1));
                [~,idxUvpDepth] = min(abs(uvpDepths - modelledTargetDepths));

                % Extract data
                modTotPnumMonthlyDh(iDh,iMonth,iLoc) = modTotPnumMonthlyProfile(idxUvpDepth,iMonth,iLoc);
                uvpTotPnumMonthlyDh(iDh,iMonth,iLoc) = uvpTotPnumMonthlyTarget(idxUvpDepth,iMonth,iLoc);
            end
        end
    end

    % Define output arrays
    if iTest == 1
        nTracers = size(modFluxMonthlyDh,5);
        arrayModFluxes = NaN(nTargetDepths,12,config.nLocs,nTracers,length(testrunDir));
        arrayModTotPnum = NaN(nTargetDepths,12,config.nLocs,length(testrunDir));
        arrayObsFluxes = NaN(nTargetDepths,12,config.nLocs,nTracers);
        arrayObsFluxesErr = NaN(nTargetDepths,12,config.nLocs,nTracers);
        arrayObsTotPnum = NaN(nTargetDepths,12,config.nLocs);
    end
    
    % Store in output array
    arrayModFluxes(:,:,:,:,iTest) = squeeze(modFluxMonthlyDh(:,:,:,1,:));
    arrayModTotPnum(:,:,:,iTest) = modTotPnumMonthlyDh;
    for iTracer = 1:nTracers
        arrayObsFluxes(:,:,:,iTracer) = obsDataMonthly{iTracer}(:,:,:);
        arrayObsFluxesErr(:,:,:,iTracer) = obsErrorMonthly{iTracer}(:,:,:);
    end
    arrayObsTotPnum = uvpTotPnumMonthlyDh;

end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - PREPARING ARRAYS TO PLOT AND CALCULATE SENSITIVITY TO 
% PARAMETER CHANGE 
% -------------------------------------------------------------------------

nParams = length(paramMetadata);
nTargetVars = 8; % POC flux zeu, POC flux zmeso, PIC flux zeu, PIC flux zmeso, bSi flux zeu, bSi flux zmeso, pnum zeu, pnum zmeso

% Dimensions: 12 months x 6 locations x (1 + nParams*2) tests + standard x 8 variables
modDataForPlotting = NaN(12,config.nLocs,(1+nParams*2),nTargetVars); 
modDataForPlotting(:,:,:,1) = squeeze(arrayModFluxes(1,:,:,1,:)); % 1 = POC at export
modDataForPlotting(:,:,:,2) = squeeze(arrayModFluxes(3,:,:,1,:)); % 2 = POC at sequestration
modDataForPlotting(:,:,:,3) = squeeze(arrayModFluxes(1,:,:,2,:)); % 3 = PIC at export
modDataForPlotting(:,:,:,4) = squeeze(arrayModFluxes(3,:,:,2,:)); % 4 = PIC at sequestration
modDataForPlotting(:,:,:,5) = squeeze(arrayModFluxes(1,:,:,3,:)); % 5 = bSi at export
modDataForPlotting(:,:,:,6) = squeeze(arrayModFluxes(3,:,:,3,:)); % 6 = bSi at sequestration
modDataForPlotting(:,:,:,7) = squeeze(arrayModTotPnum(1,:,:,:));  % 7 = Particle num at export
modDataForPlotting(:,:,:,8) = squeeze(arrayModTotPnum(3,:,:,:));  % 8 = Particle num at sequestration

% Dimensions: 12 months x 6 locations x 8 variables
obsDataForPlotting = NaN(12,config.nLocs,nTargetVars); 
obsDataForPlotting(:,:,1) = squeeze(arrayObsFluxes(1,:,:,1)); % 1 = POC at export
obsDataForPlotting(:,:,2) = squeeze(arrayObsFluxes(3,:,:,1)); % 2 = POC at sequestration
obsDataForPlotting(:,:,3) = squeeze(arrayObsFluxes(1,:,:,2)); % 3 = PIC at export
obsDataForPlotting(:,:,4) = squeeze(arrayObsFluxes(3,:,:,2)); % 4 = PIC at sequestration
obsDataForPlotting(:,:,5) = squeeze(arrayObsFluxes(1,:,:,3)); % 5 = bSi at export
obsDataForPlotting(:,:,6) = squeeze(arrayObsFluxes(3,:,:,3)); % 6 = bSi at sequestration
obsDataForPlotting(:,:,7) = squeeze(arrayObsTotPnum(1,:,:));  % 7 = Particle num at export
obsDataForPlotting(:,:,8) = squeeze(arrayObsTotPnum(3,:,:));  % 8 = Particle num at sequestration

% Dimensions: 12 months x 6 locations x 8 variables
obsErrorForPlotting = NaN(12,config.nLocs,nTargetVars); 
obsErrorForPlotting(:,:,1) = squeeze(arrayObsFluxesErr(1,:,:,1)); % 1 = POC at export
obsErrorForPlotting(:,:,2) = squeeze(arrayObsFluxesErr(3,:,:,1)); % 2 = POC at sequestration
obsErrorForPlotting(:,:,3) = squeeze(arrayObsFluxesErr(1,:,:,2)); % 3 = PIC at export
obsErrorForPlotting(:,:,4) = squeeze(arrayObsFluxesErr(3,:,:,2)); % 4 = PIC at sequestration
obsErrorForPlotting(:,:,5) = squeeze(arrayObsFluxesErr(1,:,:,3)); % 5 = bSi at export
obsErrorForPlotting(:,:,6) = squeeze(arrayObsFluxesErr(3,:,:,3)); % 6 = bSi at sequestration
obsErrorForPlotting(:,:,7) = zeros(12,config.nLocs);              % 7 = Particle num at export
obsErrorForPlotting(:,:,8) = zeros(12,config.nLocs);              % 8 = Particle num at sequestration

% Calculate effect of parameter change. Evaluate two metrics: log10 ratio
% and gradient
log10ratioSensitivityMonthly = NaN(nParams,12,config.nLocs,nTargetVars); 
log10ratioSensitivityAnnual = NaN(nParams,config.nLocs,nTargetVars);
gradientDiffMonthly = NaN(nParams,12,config.nLocs,nTargetVars); 
gradientDiffAnnual = NaN(nParams,config.nLocs,nTargetVars); 

for iVar = 1:nTargetVars
    for iLoc = 1:config.nLocs
        for iParam = 1:nParams
            paramDiff = cell2mat(paramMetadata(iParam,4)) - cell2mat(paramMetadata(iParam,3));
            idxTestLow = 2 + (iParam - 1) * 2;
            idxTestUpp = idxTestLow + 1;
            for iMonth = 1:12
                valVarAtMaxParam = modDataForPlotting(iMonth,iLoc,idxTestUpp,iVar);
                valVarAtMinParam = modDataForPlotting(iMonth,iLoc,idxTestLow,iVar);
                if valVarAtMinParam > 0 && valVarAtMaxParam > 0
                    log10ratioSensitivityMonthly(iParam,iMonth,iLoc,iVar) =...
                        log10(valVarAtMaxParam/valVarAtMinParam);
                else
                    log10ratioSensitivityMonthly(iParam,iMonth,iLoc,iVar) = NaN;
                end
                changeInTheVar = abs(valVarAtMaxParam-valVarAtMinParam);
	            gradientDiffMonthly(iParam,iMonth,iLoc,iVar) = changeInTheVar./paramDiff;
            end
            % Average over months to get annual values
            log10ratioSensitivityAnnual(iParam,iLoc,iVar) = nanmean(log10ratioSensitivityMonthly(iParam,:,iLoc,iVar));  
            gradientDiffAnnual(iParam,iLoc,iVar) = nanmean(gradientDiffMonthly(iParam,:,iLoc,iVar));
        end
    end
end

% Convert the log10-based sensitivity values into percentages
ratioMonthly = 10.^log10ratioSensitivityMonthly; % from log10 back to ratio (e.g., from 0.3 to 2)
percentageChangeMonthly = (ratioMonthly - 1) .* 100; % convert ratio to percentage change (e.g., from 2 to +100%)
ratioAnnual = 10.^log10ratioSensitivityAnnual; 
percentageChangeAnnual = (ratioAnnual - 1) .* 100;

% Take the top 10 parameters for each variable for the metric 
% log10ratioSensitivityAnnual, and only consider first 5 locations (i.e.,
% exclude anomalous HAUSGARTEN)
firstFiveLocs = squeeze(percentageChangeAnnual(:,1:5,:));
paramLabels = {paramMetadata{:,1}}; % cell array of parameter names

% First, take mean of absolute % change across locations for each parameter 
% and variable
meanAbsPctChangePerParamVar = squeeze( nanmean(abs(firstFiveLocs), 2) ); % nParams x nVars
top10ParamsPerVar = zeros(10,nTargetVars); % store indices of top 10 params
top10ValuesPerVar = zeros(10,nTargetVars);
top10ParamNamesPerVar = cell(10,nTargetVars);
for iVar = 1:nTargetVars
    [sortedVals, sortedIdx] = sort(meanAbsPctChangePerParamVar(:,iVar), 'descend');
    top10ParamsPerVar(:,iVar) = sortedIdx(1:10);
    top10ValuesPerVar(:,iVar) = sortedVals(1:10);
    top10ParamNamesPerVar(:,iVar) = paramLabels(sortedIdx(1:10))';
end

% Second, take mean of absolute % change overall
meanAbsPctChangePerParam = squeeze( nanmean(abs(firstFiveLocs), [2, 3]) ); % nParams
[sortedValues, sortedIndices] = sort(meanAbsPctChangePerParam, 'descend');
top10ParamsOverall = sortedIndices(1:10);
top10ValuesOverall = sortedValues(1:10);
top10ParamNamesOverall = paramLabels(sortedIndices(1:10))';

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - PLOT HEATMAPS
% -------------------------------------------------------------------------

paramSelect = {...
    '$$\theta_{\mathrm{bSi:C,diat}}$$', ...
    '$$\theta_{\mathrm{calc:C,cocco}}^{\mathrm{max}}$$', ...
    '$$f_{\mathrm{TEP}}$$', ...
    '$$R_{\mathrm{z}}$$', ...
    '$$\theta_{\mathrm{agg:zoo}}$$', ...
    '$$f_{\mathrm{OM}}^{\mathrm{ing,surf}}$$', ...
    '$$\lambda_{\mathrm{m,C_{POC}}}^{\mathrm{max,0}}$$', ...
    '$$\lambda_{\mathrm{m,C_{TEP}}}^{\mathrm{max,0}}$$', ...
    '$$\omega_{\mathrm{C_{TEP}}}$$', ...
    '$$\kappa_{\mathrm{bSi}}^{\mathrm{max,0}}$$' ...
}; 

varNames = {... % as chosen in section above
    'POC flux at export',...
    'POC flux at sequestration',...
    'PIC flux at export',...
    'PIC flux at sequestration',...
    'bSi flux at export',...
    'bSi flux at sequestration',...
    'Particle num. at export',...
    'Particle num. at sequestration',...
};

varNamesChosenOrder = {... % first export vars, second mesopelagic vars
    'POC flux at export',...
    'PIC flux at export',...
    'bSi flux at export',...
    'Particle num. at export',...
    'POC flux at sequestration',...
    'PIC flux at sequestration',...
    'bSi flux at sequestration',...
    'Particle num. at sequestration',...
};

% Extract original parameter names from metadata
originalParamNames = paramMetadata(:,1);

% Match selected parameters to their indices in the original list
selectedParamIdxs = zeros(numel(paramSelect),1);
for i = 1:numel(paramSelect)
    idx = find(strcmp(originalParamNames, paramSelect(i)));
    if isempty(idx)
        error('Parameter "%s" not found in paramMetadata.', paramSelect(i));
    end
    selectedParamIdxs(i) = idx;
end

% Reorder variable indices according to desired plot order
[~, reorderVarNameIdx] = ismember(varNamesChosenOrder, varNames);

% Define subplot titles (used in colour bar labels)
cbLabels = {
    {'Percentage change in', 'POC flux (%)'},...
    {'Percentage change in', 'PIC flux (%)'},...
    {'Percentage change in', 'bSi flux (%)'},...
    {'Percentage change in', 'particle number concentration (%)'},...
    {'Percentage change in', 'POC flux (%)'},...
    {'Percentage change in', 'PIC flux (%)'},...
    {'Percentage change in', 'bSi flux (%)'},...
    {'Percentage change in', 'particle', 'number concentration (%)'}...
};

% Figure configuration
layout = struct(...
    'leftMargin', 0.08, ...
    'rightMargin', 0.05, ...
    'topMargin', 0.125, ...
    'bottomMargin', 0.12, ...
    'gapX', 0.025, ...
    'gapY', 0.03, ...
    'nCols', 4, ...
    'nRows', 2 ...
);

width = (1 - layout.leftMargin - layout.rightMargin - (layout.nCols - 1) * layout.gapX) / layout.nCols;
height = (1 - layout.topMargin - layout.bottomMargin - (layout.nRows - 1) * layout.gapY) / layout.nRows;
axHandles = gobjects(layout.nRows*layout.nCols,1);
axPositions = zeros(layout.nRows*layout.nCols,4);

figure('Units', 'normalized', 'Position', [0.05, 0.1, 0.55, 0.45])

for iSubplot = 1:nTargetVars

    row = floor((iSubplot-1) / layout.nCols);
    col = mod((iSubplot-1), layout.nCols); 
    left = layout.leftMargin + col * (width + layout.gapX);
    bottom = 1 - layout.topMargin - (row + 1) * height - row * layout.gapY;
    ax = axes('Position', [left bottom width height]);
    axHandles(iSubplot) = ax;
    axPositions(iSubplot,:) = [left bottom width height];
    
    data = percentageChangeAnnual(selectedParamIdxs, reorderModLocIdx, reorderVarNameIdx(iSubplot));
    h = imagesc(ax, data);

    % Set transparency: fully opaque (1) for real numbers, transparent (0) for NaNs
    set(h, 'AlphaData', ~isnan(data));  % NaNs get 0 (transparent), others get 1

    if ismember(iSubplot, [1, 5])
        caxis(ax, [-100 100]);
    elseif ismember(iSubplot, [2, 6])
        caxis(ax, [-500 500]); 
    elseif ismember(iSubplot, [3, 7])
        caxis(ax, [-500 500]); 
    elseif ismember(iSubplot, [4, 8])
        caxis(ax, [-500 500]);
    end

    % Define which colormap to use for each subplot
    if ismember(iSubplot, [1, 5])
        currentCmap = cmocean('balance');
    elseif ismember(iSubplot, [2, 3, 6, 7])
        currentCmap = cmocean('delta');
    elseif ismember(iSubplot, [4, 8])
        currentCmap = cmocean('curl');
    end

    colormap(ax,currentCmap);

    % Only add colorbar above selected subplots
    if ismember(iSubplot, [1, 2, 3, 4])

        % Get axes position
        axPos = get(ax, 'Position');
        cb = colorbar(ax, 'Location', 'north');
        cb.XAxisLocation = 'top';  % this moves ticks to the top, so they face upwards
    
        % Match colorbar to the same color range and map
        colormap(ax,currentCmap);
    
        % Set label 
        cb.Label.String = cbLabels{iSubplot};
        cb.Label.Interpreter = 'none';  % so \n works as newline
        cb.Label.FontSize = 13;
        cb.Label.FontWeight = 'bold';
        cb.Label.Position(2) = cb.Label.Position(2) + 1.3;  % shift upward

        % Adjust position of the colorbar 
        cbPos = get(cb, 'Position');
        cbPos(1) = axPos(1);         % align left edge with axes left
        cbPos(3) = axPos(3);         % set width equal to axes width
        cbPos(2) = cbPos(2) + 0.060; % move colorbar slightly upwards (adjust as needed)
        cbPos(4) = 0.015;
        set(cb, 'Position', cbPos);

    end

    axis tight

    % Set y-ticks with latex labels
    set(ax, 'YTick', 1:length(paramSelect), ...
        'YTickLabel', paramSelect, ...
        'TickLabelInterpreter', 'latex')

    % Set x-ticks labels
    set(ax, 'XTick', 1:length(desiredModLocationOrder), ...
        'XTickLabel', desiredModLocationOrder)

    % For first row, hide x tick labels
    if ismember(iSubplot, [1, 2, 3, 4])
        ax.XTickLabel = [];
    else
        ax.XTickLabel = desiredModLocationOrder; % show for bottom row
        xtickangle(ax, 45)
    end

    % Hide Y labels except for 1st and 5th, and 4th and 8th
    if ismember(iSubplot, [4, 8]) % show y-tick labels on the right-hand side
        set(ax, 'YAxisLocation', 'right');
    elseif ~ismember(iSubplot, [1, 5]) % Hide y-tick labels for all other subplots except 1 and 5
        set(ax, 'YTickLabel', []);
    end

    % Add text labels
    threshold = 50;  % Set threshold for white text
    for iRow = 1:size(data,1)
        for iCol = 1:size(data,2)
            val = data(iRow,iCol);
            % Choose text color: white for dark backgrounds, black otherwise
            if abs(val) > threshold
                textColor = 'w';
            else
                textColor = 'k';
            end
            
            if ~isnan(val)
                absVal = abs(val);

                if absVal < 1
                    if val > 0
                        str = '<1';
                    elseif val < 0
                        str = '>-1';
                    else
                        str = '0';
                    end
            
                elseif absVal <= 1000
                    str = sprintf('%.0f', val);
            
                elseif absVal < 10000
                    bucket = floor(absVal / 1000) * 1000;
                    if val > 0
                        str = sprintf('>%d', bucket);
                    else
                        str = sprintf('<-%d', bucket);
                    end
            
                else
                    % For values >= 10,000, find the power of ten bucket:
                    exponent = floor(log10(absVal));
                    baseValue = 10^exponent;
            
                    % Determine the next lower power of 10 bucket to display
                    % Display "> 10^exponent" for values >= 10^exponent and < 10^(exponent+1)
                    if val > 0
                        str = sprintf('≥ 10^%d', exponent);
                    else
                        str = sprintf('≤ -10^%d', exponent);
                    end
                end
            else
                str = '';
            end

            text(iCol, iRow, str, ...
                'HorizontalAlignment', 'center', ...
                'Color', textColor, 'FontSize', 8);
        end
    end

    % Add annotations
    if any(iSubplot == [1, 5])
        labelAnnotationText = {'Export depth', 'Sequestration depth'};
        yAnnotationOffset = [0.04, 0.38];
        labelIndex = find([1, 5] == iSubplot);
        text(layout.leftMargin - 0.43, bottom + yAnnotationOffset(labelIndex), labelAnnotationText{labelIndex}, ...
            'Units', 'normalized', 'Rotation', 90, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontWeight', 'bold', 'FontSize', 12);
    end

    % Save the figure
    saveFigureInFolder(figureSubfolder,'heatmap_param_sensit')
   
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - PLOT GRADIENT MAPS
% -------------------------------------------------------------------------

load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES','STATION_TAGS')

monthLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

varNames = {...
    'POC flux at export depth (mg C m^{-2} d^{-1})',...
    'POC flux at sequestration depth (mg C m^{-2} d^{-1})',...
    'PIC flux at export depth (mg CaCO3 m^{-2} d^{-1})',...
    'PIC flux at sequestration depth (mg CaCO3 m^{-2} d^{-1})',...
    'bSi flux at export depth (mg opal m^{-2} d^{-1})',...
    'bSi flux at sequestration depth (mg opal m^{-2} d^{-1})',...
    'Particle num. at export depth (# part. L^{-1})',...
    'Particle num. at sequestration depth (# part. L^{-1})',...
};

varTags = {...
    'pocflux_zeu',...
    'pocflux_zmeso',...
    'picflux_zeu',...
    'picflux_zmeso',...
    'bsiflux_zeu',...
    'bsiflux_zmeso',...
    'pnum_zeu',...
    'pnum_zmeso',...
};

nSubplots = nParams*12; 
% paramsGroup1 = 1:12;       % first 12 parameters
% paramsGroup2 = 13:nParams; % remaining 11 parameters
% paramGroups = {1:12, 13:nParams};

for iVar = 1:nTargetVars
    for iLoc = 1:config.nLocs

        % for groupIdx = 1:2  % First figure (1:12), second figure (13:end)
        %     currentParams = paramGroups{groupIdx};
        %     nParamsGroup = numel(currentParams);
        %     nSubplots = nParamsGroup*nSubVars; 
    
        figure()
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.50 0.86],'Color','w') 
        haxis = zeros(nSubplots,1);
        iSubplot = 0;

        for iMonth = 1:12
            maxModel = max(max(squeeze(modDataForPlotting(:,iLoc,:,iVar)))); % across the nParams (tests) for all months
            minModel = min(min(squeeze(modDataForPlotting(:,iLoc,:,iVar)))); 
            maxObs = max(squeeze(obsDataForPlotting(:,iLoc,iVar) + obsErrorForPlotting(:,iLoc,iVar)));
            minObs = min(squeeze(obsDataForPlotting(:,iLoc,iVar) - obsErrorForPlotting(:,iLoc,iVar)));
            maxAll = max([maxModel, maxObs]);
            minAll = min([minModel, minObs]);

            for iParam = 1:nParams
            % for idx = 1:nParamsGroup
            %     iParam = currentParams(idx);

                idxTestLow = 2 + (iParam - 1) * 2;
                idxTestUpp = idxTestLow + 1;
                idxDef = 1;
    
                iSubplot = (iParam - 1) * 12 + iMonth; % for column-wise progression instead of default row-wise
                haxis(iSubplot) = subaxis(nParams,12,iSubplot,'Spacing',0.01,'Padding',0.0,'Margin',0.07);
                % iSubplot = (idx - 1) * nSubVars + iMonth; % for column-wise progression instead of default row-wise
                % haxis(iSubplot) = subaxis(nParamsGroup,nSubVars,iSubplot,'Spacing',0.01,'Padding',0.0,'Margin',0.07);
             
                % Variable value at the parameter limits (blue line - slope)
                xvals = [modDataForPlotting(iMonth,iLoc,idxTestLow,iVar),modDataForPlotting(iMonth,iLoc,idxTestUpp,iVar)]; % variable of interest
                % if isnan(xvals(1)) || isnan(xvals(2))
                %     disp('nan xvals')
                %     disp(iMonth)
                %     disp(iParam)
                % end
                yvals = [cell2mat(paramMetadata(iParam,3)); cell2mat(paramMetadata(iParam,4))]; % parameter limits
                if isnan(xvals(1)) || isnan(xvals(2)) % plot dotted blue line between known xval and default value
                    xvals = [modDataForPlotting(iMonth,iLoc,idxTestLow,iVar),modDataForPlotting(iMonth,iLoc,idxDef,iVar),modDataForPlotting(iMonth,iLoc,idxTestUpp,iVar)]; % variable of interest
                    yvals = [cell2mat(paramMetadata(iParam,3)),cell2mat(paramMetadata(iParam,2)),cell2mat(paramMetadata(iParam,4))]; 
                    plot(xvals,yvals,'Color',[0 0.4470 0.7410],'LineStyle',':','LineWidth',2.5); hold on
                else % plot continuous blue line
                    plot(xvals,yvals,'Color',[0 0.4470 0.7410],'LineWidth',1.5); hold on
                end
                
                % Observation
                middleYaxisRange = ((cell2mat(paramMetadata(iParam,4))-cell2mat(paramMetadata(iParam,3)))/2)+cell2mat(paramMetadata(iParam,3));
                scatter(obsDataForPlotting(iMonth,iLoc,iVar),middleYaxisRange,45,'o','MarkerEdgeColor','r'); hold on
                neg = obsDataForPlotting(iMonth,iLoc,iVar) - obsErrorForPlotting(iMonth,iLoc,iVar);
                pos = obsDataForPlotting(iMonth,iLoc,iVar) + obsErrorForPlotting(iMonth,iLoc,iVar);
                if obsErrorForPlotting(iMonth,iLoc,iVar) > 0
                    errorbar(obsDataForPlotting(iMonth,iLoc,iVar),middleYaxisRange,obsErrorForPlotting(iMonth,iLoc,iVar),...
                        'horizontal','Color','r','LineStyle','none','HandleVisibility','off'); hold on
                end
        
                % Default variable value
                scatter(modDataForPlotting(iMonth,iLoc,idxDef,iVar),cell2mat(paramMetadata(iParam,2)),20,'MarkerFaceColor','k'); hold on
                
                box on
                set(haxis(iSubplot),'FontSize',10)

                % Axis limits
                if (iVar > 6) % particle number variables, apply log
                    set(gca, 'XScale', 'log')
                    xlim([1 1e10])                                      
                    xticks([1 1000 1e6 1e9])                         
                    xticklabels({'1','10^3','10^6','10^9'})   
                    set(gca, 'XTickLabelRotation', 0) % prevent rotation
                else
                    xLo = nanmin([minAll, neg]) - 0.1 * nanmax([maxAll, pos]);
                    if xLo < 0
                        xLo = 0;
                    end
                    xHi = nanmax([maxAll, pos]) + 0.1 * nanmax([maxAll, pos]);
                    xlim([xLo, xHi])
                    setSmartXTicks(haxis(iSubplot))
                end

                ylim([cell2mat(paramMetadata(iParam,3)) cell2mat(paramMetadata(iParam,4))])
                numTicks = 3;
                Ly = get(gca,'YLim');
                set(gca,'YTick',linspace(Ly(1),Ly(2),numTicks))
                
                ytickformat('%.2f')
                if (iParam == 5 || iParam == 8 || iParam == 14 || iParam == 16)
                    ytickformat('%.1f')
                elseif (iParam == 9 || iParam == 10 || iParam == 12 || iParam == 13)
                    ytickformat('%.3f')
                elseif (iParam == 11 || iParam == 15)
                    ytickformat('%.4f')  
                end
    
                % Clear default titles and labels
                set(haxis(iSubplot), 'Title', [], 'XLabel', [], 'YLabel', []);
                set(haxis(iSubplot), 'TickLength', [0.06, 0.06]);
    
                % Title only on the first row (iParam == 1)
                if iParam == 1
                    title(monthLabel{iMonth}, 'FontWeight', 'bold');
                else
                    set(haxis(iSubplot), 'title', []);
                end
                
                % Y-label only on the first column (iMonth == 1)
                if iMonth == 1
                    ylabel(paramMetadata(iParam,1), 'FontWeight', 'bold', 'Interpreter', 'latex');
                else
                    set(haxis(iSubplot), 'YTickLabel', []);
                end
                
                % X-axis label: only one shared label on entire figure
                if (iParam == nParams && iMonth == ceil(12/2)) % center bottom
                    xlabel(varNames{iVar}, 'FontSize', 12, 'FontWeight', 'bold');
                end
                if (iParam ~= nParams)
                    set(haxis(iSubplot), 'XTickLabel', []);
                end
    
                % Super title for each location
                if (iParam == 1 && iMonth == ceil(12/2))
                    sgtitle(' ')  % reserve space so layout is not broken

                    % Custom title manually
                    annotation('textbox', [0 0.94 1 0.05], 'String', ...
                        sprintf('Effect of parameter variation on target output variable at %s', STATION_NAMES{iLoc}), ...
                        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 16);
                end

            end % iParam
        end % iMonth

        % Save the figure
        saveFigureInFolder(figureSubfolder,strcat('sensit_',varTags{iVar},'_',STATION_TAGS{iLoc}))
        % figureTag = sprintf('%s_part%d', varTags{iVar}, groupIdx);
        % saveFigureInFolder(figureSubfolder, strcat('sensit_', figureTag, '_', STATION_TAGS{iLoc}));

        % end % groupIdx

    end % iLoc  
end % iVar

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function setSmartXTicks(ax)
% setSmartXTicks sets up to 3 nicely spaced xticks and labels
% based on the order of magnitude of the current x-axis limits.
% Input: ax — axes handle (use gca if not specified)

    if nargin < 1
        ax = gca;
    end

    % Get current x-axis limits
    xlimVals = get(ax, 'XLim');
    xmin = xlimVals(1);
    xmax = xlimVals(2);

    % Determine data range
    range = xmax - xmin;
    
    if range == 0
        return;  % Avoid division by zero or degenerate axis
    end

    % Find a nice step size (1, 2, or 5 * 10^n)
    order = floor(log10(range));           % order of magnitude
    baseSteps = [1, 2, 5];
    bestStep = baseSteps(find((range ./ (baseSteps * 10^order)) <= 3, 1, 'first'));

    if isempty(bestStep)
        bestStep = 10;  % fallback
    end

    tickSpacing = bestStep * 10^order;

    % Create tick positions within limits
    tickStart = ceil(xmin / tickSpacing) * tickSpacing;
    tickEnd = floor(xmax / tickSpacing) * tickSpacing;
    ticks = tickStart : tickSpacing : tickEnd;

    % Apply ticks and labels
    set(ax, 'XTick', ticks);
    set(ax, 'XTickLabel', arrayfun(@num2str, ticks, 'UniformOutput', false));
end