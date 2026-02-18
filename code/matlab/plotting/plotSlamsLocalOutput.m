
% ======================================================================= %
%                                                                         %
% This script plots relevant SLAMS output at six time-series locations    %
% and compares it with observations.                                      %
%                                                                         %
%   Version 1.0 - Completed 12 Feb 2026                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath(fullfile('.','modelresources','external')))
addpath(genpath(fullfile('.','modelresources','internal'))) 
addpath(genpath(fullfile('.','data','processed'))) 
addpath(genpath(fullfile('.','data','raw'))) 
addpath(genpath(fullfile('.','data','interim'))) 
addpath(genpath(fullfile('.','code','matlab'))) 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Name of experiment
testrunDir = 'LOCALTS6';
figureSubfolder = 'LOCALTS6';

% Quick assessment?
quickDiagnostics = true;

% Choices to be applied
choiceTypeGridDomain = 2; % 1=global, 2=local, 3=global with zoom in factor

% Filename declarations
filenameSlamsOutput           = 'runsoutput.mat';
filenameSlamsBcpMetrics       = 'bcpmetrics.mat';
filenameSlamsRunGrid          = 'grid_run.mat';
filenameSlamsNumDepthLayers   = 'waterColNumDepthLayers.txt';
filenameSlamsInputData        = 'slamsForcing.mat';

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
fullpathInterimDataDir      = fullfile('.','data','interim');
fullpathProcessedDataDir    = fullfile('.','data','processed');
fullpathModelInputDataDir   = fullfile(fullpathTestDir,'modelinputdata');
fullpathModelRunsDir        = fullfile(fullpathTestDir,'modelruns');

% Load configuration parameters used in the model
config = loadModelConfigurationParameters(fullpathModelInputDataDir,...
    filenameSlamsRunGrid,filenameSlamsNumDepthLayers,choiceTypeGridDomain);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - FLUX DATA
% -------------------------------------------------------------------------

% Load modelled flux data and format it to pair it with observations
[modFluxMonthlyProfile,modFluxMonthlyDh,modFluxMonthlyDhDepths,...
modFluxSeasonalProfile,modFluxSeasonalDh,modFluxSeasonalDhDepths,... 
obsFluxDataMonthly,obsFluxErrorMonthly,obsFluxDataSeasonal,obsFluxErrorSeasonal,...
bicepExFluxMonthly,bicepExDepthMonthly,bicepExFluxSeasonal,bicepExDepthSeasonal] =...
    loadAndFormatModelledFluxDataAndObservations(choiceTypeGridDomain,config,...
        fullpathModelRunsDir,filenameSlamsOutput,filenameObsPocFlux,...
        filenameObsPicAndBsiFlux,filenameTimeseriesInformation,...
        filenameBicepExportFluxDunne,filenameBicepExportFluxHenson,...
        filenameBicepExportFluxLi,filenameZeu);

% Plot a comparison of profiles of flux data monthly (modelled vs observed)
if quickDiagnostics == false
isNormalised = false;
plotMonthlyModelledFluxVsObservations(isNormalised,figureSubfolder,...
    strcat(testrunDir,'_local_flux'),config,modFluxMonthlyProfile,...
    obsFluxDataMonthly,obsFluxErrorMonthly,bicepExFluxMonthly,bicepExDepthMonthly,...
    filenameTimeseriesInformation,filenameObsPocFlux,filenameObsPicAndBsiFlux)
end

% Plot mineral ratios
if quickDiagnostics == false
plotMonthlyObservedAndModelledMineralFluxRatios(figureSubfolder,...
    strcat(testrunDir,'_local_flux_ratios'),config,modFluxMonthlyDh,...
    obsFluxDataMonthly,filenameTimeseriesInformation)
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - PARTICLE NUMBER DATA
% -------------------------------------------------------------------------

% Load UVP5 metadata
load(fullfile('.','data','processed',filenameObsPnumUvp5),'diameterClassEdges','uvpDepths')

% Process data
[modPnumByScMonthlyProfile,modPnumByScSeasonalProfile,modPnumByVcMonthlyProfile,...
modPnumByVcSeasonalProfile,uvpPnumMonthlyProfile,~,uvpPnumMonthlyTargetValues,...
uvpPnumMonthlyTargetDepths,uvpPnumSeasonalTargetValues,uvpPnumSeasonalTargetDepths] =...
    loadAndFormatModelledParticleNumberDataAndObservations(choiceTypeGridDomain,config,...
        filenameTimeseriesInformation,filenameObsPnumUvp5,fullpathModelRunsDir,filenameSlamsOutput);

% Plot a comparison of profiles of total particle numbers monthly (modelled vs observed)
if quickDiagnostics == false
plotMonthlyModelledTotalPnumVsObservations(figureSubfolder,strcat(testrunDir,'_local_pnum_total'),...
    config,modPnumByScMonthlyProfile,uvpPnumMonthlyTargetValues,uvpPnumMonthlyTargetDepths,...
    filenameTimeseriesInformation)
end

% Plot heatmaps of particle abundance by size class and month for model outputs
isModelled = true;
pnumVar = modPnumByScMonthlyProfile;
pnumDepths = repmat(config.availImagSysDeployDepths,1,config.nLocs);
particleSizeClassEdges = config.particleDiameterClassBounds;
plotHeatmapOfLocalPnumBySizeClassAndMonth(isModelled,figureSubfolder,...
    strcat(testrunDir,'_local_pnum_heatmap_mod'),pnumVar,pnumDepths,...
    particleSizeClassEdges,filenameTimeseriesInformation)

% Plot heatmaps of particle abundance by size class and month for UVP5
% observations. For that, create a particle array array with one extra row 
% (1st dimension) to match size of modelled particle array (the diameter 
% classes by UVP data do not consider the 0.2–1 um size class range)
if quickDiagnostics == false
isModelled = false;
szOrig = size(uvpPnumMonthlyProfile);         % get size of uvpPnumMonthlyTargetValues (15 x 46 x 12 x 6)
szOrig(1) = szOrig(1) + 1;                    % increase the first dimension by 1
pnumVar = NaN(szOrig);
pnumVar(2:end,:,:,:) = uvpPnumMonthlyProfile; % copy original data into positions 2 to 16
pnumDepths = repmat(uvpDepths,1,config.nLocs);
particleSizeClassEdges = config.particleDiameterClassBounds;
plotHeatmapOfLocalPnumBySizeClassAndMonth(isModelled,figureSubfolder,...
    strcat(testrunDir,'_local_pnum_heatmap_obs'),pnumVar,pnumDepths,...
    particleSizeClassEdges,filenameTimeseriesInformation)
end

% Plot relative contribution of different ESD classes total particle number 
% and its relationship to POC flux
if quickDiagnostics == false
plotMonthlyRelativeContributionEsdToPocFlux(config,filenameTimeseriesInformation,...
    figureSubfolder,strcat(testrunDir,'_local_pnum_relatesd'),modFluxMonthlyProfile,...
    modPnumByScMonthlyProfile)
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - CORRELATION PLOTS: MODEL VS OBSERVATIONS
% -------------------------------------------------------------------------

plotCorrelationPocFluxObsVsModelled(figureSubfolder,...
    strcat(testrunDir,'_local_pocflux_correl'),modFluxMonthlyDh,...
    obsFluxDataMonthly,bicepExFluxMonthly,filenameTimeseriesInformation)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 5 - MODEL PERFORMANCE METRICS
% -------------------------------------------------------------------------

metricScoresOverall = calculateOverallModelPerformanceMetrics(config,...
    obsFluxDataMonthly,modFluxMonthlyDh,uvpPnumMonthlyTargetValues,...
    modPnumByScMonthlyProfile,uvpDepths,filenameTimeseriesInformation,...
    fullfile(fullpathProcessedDataDir,'overallModelPerformanceMetrics'));
metricScoresBicep = calculatePocExportFluxModelPerformanceMetrics(...
    obsFluxDataMonthly,bicepExFluxMonthly,modFluxMonthlyDh,filenameTimeseriesInformation,...
    fullfile(fullpathProcessedDataDir,'exportFluxModelPerformanceMetrics'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - SEASONAL FLUXES AND PARTICLE NUMBER DATA
% -------------------------------------------------------------------------

plotSeasonalModelledFluxAndPnumVsObservations(figureSubfolder,...
    strcat(testrunDir,'_local_seasonal'),config,modFluxSeasonalProfile,...
    obsFluxDataSeasonal,obsFluxErrorSeasonal,bicepExFluxSeasonal,bicepExDepthSeasonal,...
    modPnumByScSeasonalProfile,uvpPnumSeasonalTargetValues,uvpPnumSeasonalTargetDepths,...
    filenameTimeseriesInformation,filenameObsPocFlux,filenameObsPicAndBsiFlux,filenameObsPnumUvp5)

plotSeasonalRelativeContributionEsdToPocFlux(config,filenameTimeseriesInformation,...
    figureSubfolder,strcat(testrunDir,'_local_seasonal_pnum_relatesd'),...
    modFluxSeasonalProfile,modPnumByScSeasonalProfile,modPnumByScMonthlyProfile)

plotSeasonalRelativeContributionVsinkToPocFlux(config,filenameTimeseriesInformation,...
    figureSubfolder,strcat(testrunDir,'_local_seasonal_pnum_relatvsink'),...
    modFluxSeasonalProfile,modPnumByVcSeasonalProfile,modPnumByVcMonthlyProfile)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - PARTICLE ATTRIBUTES AND CORRELATION MATRIX WITH POC EXPORT
% FLUX
% -------------------------------------------------------------------------

% Process the data
[particleAttributesByScAnnualDhMetadata,...
 particleAttributesByScSeasonalDhMetadata,...
 particleAttributesByTypeAnnualDhMetadata,...
 fluxMonthlyDh,fluxSeasonalDh] =...
    loadAndFormatModelledParticleAttributeData(config,filenameTimeseriesInformation,...
        filenameObsPocFlux,fullpathModelRunsDir,filenameSlamsOutput);

% Snapshots of particle attributes
plotAverageParticleAttributes(particleAttributesByScAnnualDhMetadata,...
    particleAttributesByScSeasonalDhMetadata,fluxMonthlyDh,fluxSeasonalDh,...
    figureSubfolder,strcat(testrunDir,'_local_patt'))

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 8 - PARTICLE TYPES
% -------------------------------------------------------------------------

plotRelativeContributionParticleTypes(figureSubfolder,testrunDir,...
    'localparticletypes_fractions',config,filenameTimeseriesInformation,...
    filenameObsPocFlux,fullpathModelRunsDir,filenameSlamsOutput)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 9 - FLUX SINK TERMS
% -------------------------------------------------------------------------

plotSinkTerms(figureSubfolder,testrunDir,'localsinks_fractions',config,...
    fullpathModelRunsDir,filenameSlamsOutput,choiceTypeGridDomain)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 10 - ENCOUNTER TERMS
% -------------------------------------------------------------------------

plotEncounterTerms(figureSubfolder,testrunDir,'localencounters',config,...
    fullpathModelRunsDir,filenameSlamsOutput,choiceTypeGridDomain)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 11 - PFT PROBABILITY DISTRIBUTION
% -------------------------------------------------------------------------

plotPftProbabilityDistributions(figureSubfolder,testrunDir,'localpftprobs',config,...
    fullpathModelRunsDir,filenameSlamsOutput)

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 12 - INPUT FORCING
% -------------------------------------------------------------------------

% Load station-related information from observation compilation
load(fullfile(fullpathModelInputDataDir,filenameSlamsInputData),'mapStruct')
nLocs = size(mapStruct.monthlyNpp,2);

% Rearrange model locations to match observational order
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'};
desiredModLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
[~,reorderModLocIdx] = ismember(desiredModLocationOrder,currentModLocationOrder); % get reordering indices
stationTags = {'hotaloha','batsofp','eqpac','papso','osp','hausgarten'};

% Take the indices of the first step of each day in a 1095 vector (1 day
% ahs 3 time steps)
firstOfDayIdx = 1:3:1095;  

% Compute the cumulative start day for each month
daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31]; % month lengths
monthStart = cumsum([1, daysInMonth(1:end-1)]);  

% Organise data into a structure
dailyDepthDistribData = struct(...
    'Mesozoo', mapStruct.monthlyMesozoo(:,firstOfDayIdx,reorderModLocIdx), ...
    'Nit', mapStruct.monthlyNit(:,firstOfDayIdx,reorderModLocIdx), ...
    'Sil', mapStruct.monthlySil(:,firstOfDayIdx,reorderModLocIdx), ...
    'Phos', mapStruct.monthlyPhos(:,firstOfDayIdx,reorderModLocIdx), ...
    'OmegaC', mapStruct.monthlyOmegaC(:,firstOfDayIdx,reorderModLocIdx), ...
    'Oxy', mapStruct.monthlyOxy(:,firstOfDayIdx,reorderModLocIdx), ...
    'Temp', mapStruct.monthlyTemp(:,firstOfDayIdx,reorderModLocIdx), ...
    'Rho', mapStruct.monthlyRho(:,firstOfDayIdx,reorderModLocIdx), ... 
    'DynVisco', mapStruct.monthlyDynVisco(:,firstOfDayIdx,reorderModLocIdx));
dailySurfaceData = struct(...
    'Npp', mapStruct.monthlyNpp(firstOfDayIdx,reorderModLocIdx).*3600.*24.*config.molarMassCarbon.*1e3,... % mol C m-2 s-1 --> mg C m-2 d-1
    'Chla', mapStruct.monthlyChla(firstOfDayIdx,reorderModLocIdx), ...
    'Par0', mapStruct.dailyPar0daylight(:,reorderModLocIdx), ...
    'Mld', mapStruct.monthlyMld(firstOfDayIdx,reorderModLocIdx),...
    'Dust', mapStruct.monthlyDust(firstOfDayIdx,reorderModLocIdx).*3600.*24.*1e3); % g clay m-2 s-1 --> mg clay m-2 d-1

% Datasets metadata
dailyDepthDistribDataMetadata = {...  
    'Mesozoo', sprintf('Mesozoo. (mg C m^{-3})'),             0,   6,; 
    'Nit',     sprintf('Nitrate (mmol m^{-3})'),              0,   35;        
    'Sil',     sprintf('Silicate (mmol m^{-3})'),             0,   70;  
    'Phos',    sprintf('Phosphate (mmol m^{-3})'),            0,  2.6;                          
    'OmegaC',  sprintf('Omega calcite'),                      2,    6;  
    'Oxy',     sprintf('Oxygen (mL L^{-1})'),                 2,    8; 
    'Temp',    sprintf('Temperature (ºC)'),                  -2,   30; 
    'Rho',     sprintf('Density (g cm^{-3})'),             1.02, 1.03;
    'DynVisco',sprintf('Dyn. visco. (g cm^{-1} s^{-1})'), 0.008, 0.02;        
};
dailySurfaceDataMetadata = {... 
    'Npp',     sprintf('NPP (mg C m^{-2} s^{-1})'),       200, 1050;
    'Chla',    sprintf('Chl a (mg m^{-3})'),                0, 1.05; 
    'Par0',    sprintf('PAR0 (W m^{-2})'),                  0,  180;
    'Mld',     sprintf('MLD (m)'),                          0,  270; 
    'Dust',    sprintf('Clay (mg m^{-2} d^{-1})'),       1e-2,  100;          
};

plotLocalForcingDataDepthDistributed(dailyDepthDistribData,dailyDepthDistribDataMetadata,...
    config,monthStart,desiredModLocationOrder,stationTags,figureSubfolder,testrunDir)

plotLocalForcingDataSurfaceDistribution(dailySurfaceData,dailySurfaceDataMetadata,...
    config,monthStart,desiredModLocationOrder,stationTags,figureSubfolder,testrunDir)

% % =========================================================================
% %%
% % -------------------------------------------------------------------------
% % SECTION 13 - INSTANTANEOUS SNAPSHOTS OF PARTICLE ATTRIBUTES
% % -------------------------------------------------------------------------
% 
% %% SECTION 1: Load model output data
% 
% load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')
% 
% % 20,000 x 26 x nTimeSteps in last year x nLocs
% tseptSnapshots = output.particle.snapshots;
% nSnapshots = size(tseptSnapshots, 3);
% 
% % Rearrange model locations
% currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'};
% desiredModLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
% [~,reorderModLocIdx] = ismember(desiredModLocationOrder,currentModLocationOrder); % get reordering indices
% 
% % Reorder to match desired locations
% tseptSnapshots = tseptSnapshots(:,:,:,reorderModLocIdx);
% 
% %% Start simple, with one location
% iLoc = 1;
% 
% radius   = squeeze(output.particle.snapshots(:,config.cluster.radius,:,iLoc));    % [clusters × time]
% velocity = squeeze(output.particle.snapshots(:,config.cluster.velocity,:,iLoc));
% depth    = squeeze(output.particle.snapshots(:,config.cluster.depth,:,iLoc));
% phase    = squeeze(output.particle.snapshots(:,config.cluster.phase,:,iLoc));
% nPxC     = squeeze(output.particle.snapshots(:,config.cluster.nPxC,:,iLoc));
% nPxC = double(nPxC);
% 
% % Mask --> create one long list of particles over the year, properly weighted
% mask = (phase == 1) & (nPxC > 0);
% radius   = radius(mask);
% velocity = velocity(mask);
% depth    = depth(mask);
% weights  = nPxC(mask);
% 
% %% Particle size distribution
% diam = 2*radius;  % um
% 
% mask = isfinite(diam) & isfinite(weights) & weights > 0;
% diam    = diam(mask);
% weights = weights(mask);
% 
% edges = config.particleDiameterClassBounds(:)';
% edges = unique(edges, 'stable');
% 
% binIdx = discretize(diam, edges);
% 
% nBins = numel(edges) - 1;
% 
% countSize = accumarray( ...
%     binIdx(~isnan(binIdx)), ...
%     weights(~isnan(binIdx)), ...
%     [nBins 1], ...
%     @sum, ...
%     0);
% 
% if sum(countSize) > 0
%     pdfSize = countSize / sum(countSize);
% else
%     pdfSize = zeros(size(countSize));
% end
% 
% 
% figure
% stairs(edges(1:end-1), pdfSize, 'LineWidth', 2)
% set(gca,'XScale','log')
% xlabel('Particle diameter (\mum)')
% ylabel('Relative abundance')
% title('Annual particle size distribution')
% 
% 
% %% Particle velocity distribution
% [countVel, edgesVel] = histcounts(velocity,config.particleVelocityClassBounds,'Weights',weights);
% pdfVel = countVel / sum(countVel);
% 
% figure
% stairs(edgesVel(1:end-1), pdfVel, 'LineWidth',2)
% xlabel('Settling velocity (m d^{-1})')
% ylabel('Relative abundance')
% title('Annual velocity distribution')
% 
% %% Loop over months and depths to build PSDs
% targetDepths = [0 50 100 200 500 1000];
% depthBuffer  = 10;                    % ±10 m buffer
% nMonths = 12;
% nDepths = numel(targetDepths);
% 
% t = 1:nSnapshots;
% dayOfYear = ceil(t/config.nTimestepsPerDay);
% timeMonth = zeros(size(dayOfYear));
% for m = 1:12
%     if m == 1
%         mask = dayOfYear <=  config.cumDays(m);
%     else
%         mask = dayOfYear > config.cumDays(m-1) & dayOfYear <= config.cumDays(m);
%     end
%     timeMonth(mask) = m;
% end
% 
% PSD_month_depth = NaN(nMonths, nDepths, numel(diamBins)-1);
% for iMonth = 1:12
% 
%     % Time mask for this month
%     tmask = (timeMonth == iMonth);
% 
%     for iZ = 1:nDepths
% 
%         z0 = targetDepths(iZ);
%         zmin = z0 - depthBuffer;
%         zmax = z0 + depthBuffer;
% 
%         % Combined mask
%         mask = ...
%             tmask & ...
%             (phase == 1) & ...
%             (depth >= zmin) & (depth <= zmax) & ...
%             (nPxC > 0);
% 
%         if ~any(mask(:))
%             continue
%         end
% 
%         diamSel   = diam(mask);
%         weights   = nPxC(mask);
% 
%         % Weighted PSD
%         counts = histcounts(diamSel, config.particleDiameterClassBounds, 'Weights', weights);
% 
%         % Normalise to relative abundance
%         PSD_month_depth(iMonth,iZ,:) = counts / sum(counts);
%     end
% end
% 
% iZ = find(targetDepths == 100);
% 
% figure; hold on
% for iMonth = 1:12
%     stairs(config.particleDiameterClassBounds(1:end-1), squeeze(PSD_month_depth(iMonth,iZ,:)), ...
%            'LineWidth', 1.5);
% end
% % set(gca,'XScale','log')
% xlabel('Particle diameter (\mum)')
% ylabel('Relative abundance')
% title('Monthly PSD at 100 m')
% legend(month(datetime(2024,1:12,1),'shortname'),'Location','eastoutside')
% 
% 
% %% Depth-resolved size distribution
% pdfSizeZ = NaN(length(depthBins),size(config.particleDiameterClassBounds));
% for iZ = 1:numel(depthBins)-1
%     zmask = depth >= depthBins(iZ) & depth < depthBins(iZ+1);
%     [cnt,~] = histcounts(diam(zmask),config.particleDiameterClassBounds,'Weights',weights(zmask));
%     pdfSizeZ(iZ,:) = cnt / sum(cnt);
% end
% 
% %% Seasonal distribution
% pdfSeason = NaN();
% for iSeason = 1:4
%     tmask = ismember(timeMonth, config.seasonIndices{iSeason});
%     diamS   = diam(:,tmask);
%     wS      = weights(:,tmask);
% 
%     diamS   = diamS(:);
%     wS      = wS(:);
% 
%     [cnt,~] = histcounts(diamS, sizeBins, 'Weights', wS);
%     pdfSeason(iSeason,:) = cnt / sum(cnt);
% end
% 
% % Joint size–velocity distribution
% [count2D, sizeEdges, velEdges] = histcounts2( ...
%     diam, velocity, config.particleDiameterClassBounds, config.particleVelocityClassBounds, ...
%     'Weights', weights);
% pdf2D = count2D / sum(count2D(:));
% 
% figure
% imagesc(sizeEdges(1:end-1), velEdges(1:end-1), log10(pdf2D'))
% axis xy
% % set(gca,'XScale','log')
% xlabel('Diameter (\mum)')
% ylabel('Velocity (m d^{-1})')
% colorbar
% title('Joint size–velocity distribution (log_{10})')

% % %% SECTION 3: Depth selection for plotting
% % 
% % POC flux compilation depths
% load(fullfile('.','data','processed',filenameObsPocFlux),'LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED')
% 
% % Rearrange observed locations depths
% [~,reorderObsLocIdx] = ismember(desiredModLocationOrder,STATION_NAMES); % get reordering indices
% LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(:,reorderObsLocIdx,:,:);
% 
% % Calculate annual average
% obsDepthsAnnualAvg = squeeze(mean(LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED,1,'omitnan'));
% obsDepthsMonthlyAvg = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED;
%
% % Calculate seasonal average
% obsDepthsSeasonalAvg = NaN(4,size(obsDepthsAnnualAvg,1),size(obsDepthsAnnualAvg,2),size(obsDepthsAnnualAvg,3));
% for iSeason = 1:4
%     obsDepthsSeasonalAvg(iSeason,:,:,:) = mean(obsDepthsMonthlyAvg(idxSeason{iSeason},:,:,:),1,'omitnan');
% end
% 
% % Replace the last observed depth in sediment traps with a depth more
% % relevant for particle number observations (last non-NaN)
% localImagedDepths = output.deployment.imagSys;
% for iLoc = 1:config.nLocs
%     idxLastObsDepth = find(~isnan(localImagedDepths(:,iLoc)), 1, 'last'); % last non-NaN
%     obsDepthsAnnualAvg(iLoc,:,4) = localImagedDepths(idxLastObsDepth,iLoc);
%     obsDepthsMonthlyAvg(:,iLoc,:,4) = localImagedDepths(idxLastObsDepth,iLoc);
%     obsDepthsSeasonalAvg(:,iLoc,:,4) = localImagedDepths(idxLastObsDepth,iLoc);
% end
% 
% %% SECTION 4: Extract modelled particle attribute data by main size class at key depth horizons
% 
% % Extract ALL attribute data and reorder
% particleSnapshots = output.particle.snapshots; % nParticles sizes x nAttributes x nLocs
% particleSnapshots = particleSnapshots(:,:,reorderModLocIdx);
% 
% attsIds = {...
%     'id',...
% 	'initType',...
% 	'phase',...
% 	'living',...
% 	'faecal',...
% 	'tstepCreat',...
% 	'molesOrgC',...
% 	'molesTepC',...
% 	'molesMineralOpal',...
% 	'molesMineralCalcite',...
% 	'molesMineralClay',...
% 	'massOrgMatter',...
% 	'massTep',...
% 	'nPxC',...
% 	'nPpxP',...
% 	'nPpxC',...
% 	'depth',...
% 	'mass',...
% 	'solidVolume',...
% 	'density',...
% 	'radius',...
% 	'radiusPp',...
% 	'porosity',...
% 	'stickiness',...
% 	'velocity',...
% };
% 
% %%
% 
% 
% % Find the attribute indices
% depthIdx = find(strcmp(attsIds,'depth'));
% wIdx     = find(strcmp(attsIds,'velocity'));
% 
% depth_targets = [100, 300, 500, 1000];  % horizons
% tol = 10;   % +/- tolerance window in meters
% 
% for iLoc = 1%:config.nLocs
%     % Example: extract depth and sinking velocity at location 1
%     depths = particleSnapshots(:,depthIdx,iLoc);
%     vels   = particleSnapshots(:,wIdx,iLoc);
% 
%     % Depth horizons of interest
%     depthHorizons = [50, 150, 500, 1000];
%     tol = 10; % tolerance around each horizon (m)
% 
%     figure;
%     for i = 1:numel(depthHorizons)
%         targetDepth = depthHorizons(i);
% 
%         % Select particles within tolerance of the depth horizon
%         mask = abs(depths - targetDepth) <= tol & ~isnan(vels);
%         velsSubset = vels(mask);
% 
%         % Plot histogram in subplot
%         subplot(2,2,i);
%         histogram(velsSubset, 'BinWidth', 1); % adjust bin width if needed
%         xlabel('Sinking velocity (m d^{-1})');
%         ylabel('Particle count');
%         title(sprintf('Depth: %d m (±%d m)', targetDepth, tol));
%     end
%     sgtitle(sprintf('Velocity distributions at location %d', iLoc));
% end
% 
% %%
% 
% % Extract particle attribute data at depths of interest
% nTargetDepths = size(obsDepthsAnnualAvg,3);
% nImagedDepths = size(particleAvgAttsAnnualProfile,2);
% nPartAtts = size(particleAvgAttsAnnualProfile,3);
% nPartSizes = size(particleAvgAttsAnnualProfile,1);
% particleAvgAttsAnnualDh = NaN(nPartSizes,nTargetDepths,nPartAtts,config.nLocs);
% particleAvgAttsSeasonalDh = NaN(nPartSizes,nTargetDepths,nPartAtts,4,config.nLocs);
% 
% validDepthMaskAnnual = false(nPartSizes,nImagedDepths,nPartAtts,config.nLocs,nTargetDepths);
% validDepthMaskSeasonal = false(nPartSizes,nImagedDepths,nPartAtts,4,config.nLocs,nTargetDepths);
% for iDh = 1:nTargetDepths
%     for iLoc = 1:config.nLocs
%         depthRangeTop = obsDepthsAnnualAvg(iLoc,1,iDh);    % upper boundary depth
%         depthRangeBottom = obsDepthsAnnualAvg(iLoc,2,iDh); % lower boundary depth
% 
%         % Only proceed if both are numbers rather than NaN
%         if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
%             % Find the closest depth greater than or equal to depthRangeTop
%             [~,idxTop] = min(abs(localImagedDepths(:,iLoc) - depthRangeTop));
% 
%             % Find the closest depth less than or equal to depthRangeBottom
%             [~,idxBottom] = min(abs(localImagedDepths(:,iLoc) - depthRangeBottom));
%             if idxBottom > nImagedDepths
%                 idxBottom = nImagedDepths;
%             end
% 
%             % Extract valid modelled depth range
%             validDepthMaskAnnual(:,(idxTop:idxBottom),:,iLoc,iDh) = true;
%             validDepthMaskSeasonal(:,(idxTop:idxBottom),:,:,iLoc,iDh) = true;
%         end
%     end % iLoc
% end % iDh
% 
% % Extract particle attribute data
% for iDh = 1:nTargetDepths 
%     % Annual
%     currentMaskAnnual = squeeze(validDepthMaskAnnual(:,:,:,:,iDh));
%     partAttValuesAnnual = particleAvgAttsAnnualProfile; 
%     partAttValuesAnnual(~currentMaskAnnual) = NaN;
%     particleAvgAttsAnnualDh(:,iDh,:,:) = squeeze(mean(partAttValuesAnnual,2,'omitnan'));
%     % Seasonal
%     currentMaskSeasonal = squeeze(validDepthMaskSeasonal(:,:,:,:,:,iDh));
%     partAttValuesSeasonal = particleAvgAttsSeasonalProfile; 
%     partAttValuesSeasonal(~currentMaskSeasonal) = NaN;
%     particleAvgAttsSeasonalDh(:,iDh,:,:,:) = squeeze(mean(partAttValuesSeasonal,2,'omitnan'));
% end
% 
% %% SECTION 5: Extract modelled particle attribute data by particle type at key depth horizons
% 
% % Extract ALL attribute data and reorder
% particleAvgAttInTyProfile = output.particle.avgAttInTyAnnual; % 8 types x nDepths x nAtts x nLocs
% particleAvgAttInTyProfile = particleAvgAttInTyProfile(:,:,:,reorderModLocIdx);
% 
% % Extract particle attribute data at depths of interest
% nPartTypes = size(particleAvgAttInTyProfile,1);
% nPartAtts = size(particleAvgAttInTyProfile,3);
% particleAvgAttInTyAnnualDh = NaN(nPartTypes,nTargetDepths,nPartAtts,config.nLocs);
% particleAvgAttInTyAnnualDhDepths = NaN(size(particleAvgAttInTyAnnualDh));
% 
% validDepthMaskAnnual = false(nPartTypes,nImagedDepths,nPartAtts,config.nLocs,nTargetDepths);
% for iDh = 1:nTargetDepths
%     for iLoc = 1:config.nLocs
%         depthRangeTop = obsDepthsAnnualAvg(iLoc,1,iDh);    % upper boundary depth
%         depthRangeBottom = obsDepthsAnnualAvg(iLoc,2,iDh); % lower boundary depth
% 
%         % Only proceed if both are numbers rather than NaN
%         if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
%             % Find the closest depth greater than or equal to depthRangeTop
%             [~,idxTop] = min(abs(localImagedDepths(:,iLoc) - depthRangeTop));
% 
%             % Find the closest depth less than or equal to depthRangeBottom
%             [~,idxBottom] = min(abs(localImagedDepths(:,iLoc) - depthRangeBottom));
%             if idxBottom > nImagedDepths
%                 idxBottom = nImagedDepths;
%             end
% 
%             % Extract valid modelled depth range
%             validDepthMaskAnnual(:,(idxTop:idxBottom),:,iLoc,iDh) = true;
% 
%             % Calculate the average monthly unique modelled depth
%             particleAvgAttInTyAnnualDhDepths(:,iDh,:,iLoc) =...
%                 (localImagedDepths(idxTop,iLoc)+localImagedDepths(idxBottom,iLoc))/2;
%         end
%     end % iLoc
% end % iDh
% 
% % Extract particle attribute data
% for iDh = 1:nTargetDepths
%     currentMaskAnnual = squeeze(validDepthMaskAnnual(:,:,:,:,iDh));
%     partAttValuesAnnual = particleAvgAttInTyProfile(:,:,:,:); 
%     partAttValuesAnnual(~currentMaskAnnual) = NaN;
%     particleAvgAttInTyAnnualDh(:,iDh,:,:) = squeeze(mean(partAttValuesAnnual,2,'omitnan'));
% end
% 
% %% SECTION 6: Extract modelled flux components at key depth horizons
% 
% % Extract flux components
% tracers = struct(...
%     'POC', squeeze(output.flux.monthly(1,:,:,:)) + squeeze(output.flux.monthly(2,:,:,:)),...
%     'PIC', squeeze(output.flux.monthly(3,:,:,:)),...
%     'BSi', squeeze(output.flux.monthly(4,:,:,:))...
% );
% tracerNames = fieldnames(tracers);
% 
% % Assemble modelled flux array (nDepths x 12 x nLocs x 3 tracers)
% fluxMonthlyProfile = NaN([size(tracers.POC),numel(tracerNames)]); 
% for iTracer = 1:numel(tracerNames)
%     fluxMonthlyProfile(:,:,:,iTracer) = tracers.(tracerNames{iTracer}); 
% end 
% fluxMonthlyProfile = fluxMonthlyProfile(:,:,reorderModLocIdx,:);
% 
% % Extract flux data at depths of interest
% localSedTrapDepths = output.deployment.sedTrap;
% nTargetDepths = size(obsDepthsMonthlyAvg,4);
% nSedTrapDepths = size(fluxMonthlyProfile,1);
% fluxMonthlyDh = NaN(nTargetDepths,12,config.nLocs,numel(tracerNames));
% fluxMonthlyDhDepths = NaN(size(fluxMonthlyDh));
% 
% validDepthMaskMonthly = false(nSedTrapDepths,12,config.nLocs,numel(tracerNames),nTargetDepths);
% for iDh = 1:nTargetDepths
%     for iMonth = 1:12
%         for iLoc = 1:config.nLocs
%             depthRangeTop = obsDepthsMonthlyAvg(iMonth,iLoc,1,iDh);    % upper boundary depth
%             depthRangeBottom = obsDepthsMonthlyAvg(iMonth,iLoc,2,iDh); % lower boundary depth
% 
%             % Only proceed if both are numbers rather than NaN
%             if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
%                 % Find the closest depth greater than or equal to depthRangeTop
%                 [~,idxTop] = min(abs(localSedTrapDepths(:,iLoc) - depthRangeTop));
% 
%                 % Find the closest depth less than or equal to depthRangeBottom
%                 [~,idxBottom] = min(abs(localSedTrapDepths(:,iLoc) - depthRangeBottom));
%                 if idxBottom > nSedTrapDepths
%                     idxBottom = nSedTrapDepths;
%                 end
% 
%                 % Extract valid modelled depth range
%                 validDepthMaskMonthly((idxTop:idxBottom),iMonth,iLoc,:,iDh) = true;
% 
%                 % Calculate the average monthly unique modelled depth
%                 fluxMonthlyDhDepths(iDh,iMonth,iLoc,:) =...
%                     (localSedTrapDepths(idxTop,iLoc)+localSedTrapDepths(idxBottom,iLoc))/2;
%             end
%         end % iLoc
%     end % iMonth
% end % iDh
% 
% % Extract flux data
% for iDh = 1:nTargetDepths
%     currentMask = squeeze(validDepthMaskMonthly(:,:,:,:,iDh));
%     fluxValues = fluxMonthlyProfile; 
%     fluxValues(~currentMask) = NaN;
%     fluxMonthlyDh(iDh,:,:,:) = squeeze(mean(fluxValues,1,'omitnan'));
% end
% 
% % Calculate seasonal average
% fluxSeasonalDh = NaN(size(fluxMonthlyDh,1),4,config.nLocs,numel(tracerNames));
% for iSeason = 1:4
%     fluxSeasonalDh(:,iSeason,:,:) = mean(fluxMonthlyDh(:,idxSeason{iSeason},:,:),2,'omitnan');
% end
% 
% %% SECTION 7: Assemble particle attribute into a structure for easier management of the array
% 
% % Got rid of attribute #9, depth
% particleAttributesByScAnnualDhMetadata = {... % nLocs x 2 sizes x 4 depths
%     'poc',   '$${\mathrm{POC}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,10,:)), [3,1,2]), 'POC content (pmol)';
%     'tepc',  '$${\mathrm{TEP-C}}$$',            permute(squeeze(particleAvgAttsAnnualDh(2:3,:,11,:)), [3,1,2]), 'TEP-C content (pmol)';
%     'alpha', '$$\alpha_{\mathrm{p}}$$',         permute(squeeze(particleAvgAttsAnnualDh(2:3,:,4,:)), [3,1,2]),  'Stickiness';
%     'calc',  '$${\mathrm{CaCO_{3}}}$$',         permute(squeeze(particleAvgAttsAnnualDh(2:3,:,13,:)), [3,1,2]), 'CaCO3 content (pmol)';
%     'bsi',   '$${\mathrm{bSi}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,12,:)), [3,1,2]), 'Opal content (pmol)';
%     'clay',  '$${\mathrm{Clay}}$$',             permute(squeeze(particleAvgAttsAnnualDh(2:3,:,14,:)), [3,1,2]), 'Clay content (pmol)';    
%     'vol',   '$$V_{\mathrm{p}}^{\mathrm{S}}$$', permute(squeeze(particleAvgAttsAnnualDh(2:3,:,6,:)), [3,1,2]),  'Solid volume (\mu^{3})';
%     'd3',    '$$D_{\mathrm{3}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,7,:)), [3,1,2]),  'Fractal dimension';
%     'esd',   '$${\mathrm{ESD}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,15,:)), [3,1,2]), 'ESD (\mum)';
%     'por',   '$$P_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,5,:)), [3,1,2]),  'Porosity';
%     'rho',   '$$\rho_{\mathrm{p}}$$',           permute(squeeze(particleAvgAttsAnnualDh(2:3,:,2,:)), [3,1,2]),  'Density (g cm^{-3})';
%     'vsink', '$$v_{\mathrm{p}}$$',              permute(squeeze(particleAvgAttsAnnualDh(2:3,:,3,:)), [3,1,2]),  'Sinking velocity (m d^{-1})';
%     'pnum',  '$$N$$',                           permute(squeeze(particleAvgAttsAnnualDh(2:3,:,1,:)), [3,1,2]),  'Particle number (# L^{-1})';
% };
% 
% % =========================================================================
% %%
% % -------------------------------------------------------------------------
% % PHYTOPLANKTON CARBON CONCENTRATION
% % -------------------------------------------------------------------------
% 
% %% SECTION 1: Load metadata
% 
% % Load station-related information from observation compilation
% load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
% nLocs = length(STATION_NAMES);
% 
% % Rearrange model locations to match observational order
% currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'}; % as in config.gridLats, config.gridLons
% desiredModLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'}; % as in config.gridLats, config.gridLons
% [~,reorderModLocIdx] = ismember(desiredModLocationOrder,currentModLocationOrder); % get reordering indices
% 
% %% SECTION 2: Load BICEP data
% 
% load(fullfile(fullpathInterimDataDir,filenameBicepPft),'pft_bicep','pft_bicep_lon','pft_bicep_lat')
% globalBicepCphytoMonthly = pft_bicep(:,:,:,2:4); % the 1st group is total phytoplankton biomass; 2, 3, 4 are Cmicro, Cnano and Cpico
% globalBicepCphytoAnnual = squeeze(mean(globalBicepCphytoMonthly,3,'omitnan')); % mg C m-3
% 
% % Extract data at our locations of interest
% localBicepCphytoAnnual = NaN(config.nLocs,3); % mg C m-3
% for iLoc = 1:config.nLocs
%     [~,iLon] = min(abs(pft_bicep_lon(:) - config.lons(iLoc)));
%     [~,iLat] = min(abs(pft_bicep_lat(:) - config.lats(iLoc)));
%     localBicepCphytoAnnual(iLoc,:) = globalBicepCphytoAnnual(iLat,iLon,:);
% end    
% localBicepCphytoAnnual = localBicepCphytoAnnual(reorderModLocIdx,:);
% 
% %% SECTION 3: Load modelled data
% 
% load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')
% tmpOutput = extractModelOutput(output);
% 
% % If loaded model output is global, extract local data
% if choiceTypeGridDomain == 1
%     tmpOutput = extractLocalOutputDataFromGlobalArray(config,nLocs,...
%         config.gridLons,config.gridLats,tmpOutput);
% end
% 
% % Assemble auxilliary information of phytoplankton mol C m-3
% phytoIdx = [
%     config.aux.idxDiatBiomass;
%     config.aux.idxFlagelBiomass;
%     config.aux.idxCoccoBiomass;
%     config.aux.idxPicoBiomass
% ];
% 
% % Extract surface (depth=1) biomass and convert units (mol C m-3 --> mg C m-3)
% modAuxAnnual = tmpOutput.annualAux; % # idxs x nDepths x nLocs
% modAuxAnnual = modAuxAnnual(:,:,reorderModLocIdx);
% modCphytoAnnual = squeeze(modAuxAnnual(phytoIdx,1,:)).*config.molarMassCarbon.*1e3; % 4 phyto types x nLocs
% 
% %% Plot bar chart
% 
% myColourScheme = brewermap(config.nLocs,'*Set1');
% 
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.40],'Color','w')
% haxis = zeros(2,1);
% 
% for iSubplot = 1:2
% 
%     haxis(iSubplot) = subaxis(2,1,iSubplot,'Spacing',0.02,'Padding',0.04,'Margin',0.08);
%     ax(iSubplot).pos = get(haxis(iSubplot),'Position');
%     set(haxis(iSubplot), 'Position', ax(iSubplot).pos)
% 
%     % Extract data
%     if iSubplot == 1 % BICEP
%         barData = localBicepCphytoAnnual'; % nLocs x 3 --> 3 x nLocs
%         tickLabels = {'Microphyto', 'Nanophyto', 'Picophyto'};
%         titleStr = 'BICEP Cphyto data';
%     else % SLAMS2.0
%         barData = modCphytoAnnual; % 4 x nLocs
%         tickLabels = {'Diatoms','Flagellates','Coccolithophores','Picophytoplankton'};
%         titleStr = 'SLAMS-2.0 output';
%     end
% 
%     % Create stacked bar plot
%     h = bar(barData);
% 
%     % Set colours for each bar group
%     for k = 1:length(h)
%         h(k).FaceColor = myColourScheme(k, :); % assign custom color to each group
%     end
%     box on
% 
%     % Add labels
%     xticks(1:length(tickLabels));
%     xticklabels(tickLabels);
%     ylabel('mg C m^{-3}');
%     title(titleStr,'FontSize',16);
% 
%     grid on
% 
% end

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function plotMonthlyObservedAndModelledMineralFluxRatios(figureSubfolderName,...
    figurePrefixName,config,modFlux,obsFlux,filenameTimeseriesInformation)

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
nLocs = length(STATION_NAMES);

% Reorder locations modelled and observed
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
currentModLocationOrder = STATION_NAMES;
[~,reorderModLocIdx] = ismember(desiredLocationOrder,currentModLocationOrder); % get reordering indices

% Arrays to plot
obsPicToPoc = NaN(size(obsFlux{1})); 
obsBSiToPoc = NaN(size(obsPicToPoc));
modPicToPoc = NaN(size(obsFlux{1})); 
modBSiToPoc = NaN(size(modPicToPoc));

% Notice units and transform
% mg C m-2 d-1 --> mol POC m-2 d-1
% mg CaCO3 m-2d-1 --> mol PIC m-2 d-1
% mg opal m-2 d-1 --> mol Si m-2 d-1
for iRatio = 1:2
    obsPocFlux = obsFlux{1}(:,:,reorderModLocIdx)./config.molarMassCarbon;
    obsPicFlux = obsFlux{2}(:,:,reorderModLocIdx)./config.molarMassCaCO3;
    obsBSiFlux = obsFlux{3}(:,:,reorderModLocIdx)./config.molarMassBiogenicSilica;
    modPocFlux = modFlux(:,:,reorderModLocIdx,1,1)./config.molarMassCarbon;
    modPicFlux = modFlux(:,:,reorderModLocIdx,1,2)./config.molarMassCaCO3;
    modBSiFlux = modFlux(:,:,reorderModLocIdx,1,3)./config.molarMassBiogenicSilica;
    if iRatio == 1
        obsPicToPoc(:,:,:) = obsPicFlux./obsPocFlux;
        modPicToPoc(:,:,:) = modPicFlux./modPocFlux;
    elseif iRatio == 2 
        obsBSiToPoc(:,:,:) = obsBSiFlux./obsPocFlux;
        modBSiToPoc(:,:,:) = modBSiFlux./modPocFlux;
    end
end

% Group data
obsData = {obsPicToPoc,obsBSiToPoc};
modData = {modPicToPoc,modBSiToPoc};

% Plot
nTargetDepths = size(obsPicToPoc,1);
flipMyColours = flipud(parula(nTargetDepths)); % flip upside down

for iRatio = 1:2
    for iType = 1:2
        if iType == 1
            valsToPlot = obsData{iRatio};
        elseif iType == 2
            valsToPlot = modData{iRatio};
        end

        figure()
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.70 0.50],'Color','w') 
        haxis = zeros(nTargetDepths,nLocs);
        iSubplot = 0;
        
        for iDh = 1:nTargetDepths
            for iLoc = 1:nLocs
                iSubplot = iSubplot + 1;
        
                haxis(iSubplot) = subaxis(nTargetDepths,nLocs,iSubplot,'Spacing',0.01,'Padding',0.01,'Margin', 0.09);
                ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        
                if (iSubplot >= 1 && iSubplot <= 6)
                    ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.04;
                elseif (iSubplot >= 6 && iSubplot <= 12)
                    ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.01; 
                elseif (iSubplot >= 12 && iSubplot <= 18)
                    ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.02;
                elseif (iSubplot >= 18)
                    ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.05;
                end
                set(haxis(iSubplot),'Position',ax(iSubplot).pos)
        
                hbar = bar(haxis(iSubplot),(1:12),squeeze(valsToPlot(iDh,:,iLoc)),...
                    'BarWidth',0.75,'FaceColor','flat');
                hbar.CData(:,:) = repmat(flipMyColours(iDh,:),[12 1]);
                hold off
                
                % Base of euphotic
                if (iSubplot >= 1 && iSubplot <= 6)
                    if iType == 1 && iRatio == 1
                        ylim([0 1.6])
                        yTickValues = 0:0.5:1.5;
                    elseif iType == 2 && iRatio == 1
                        ylim([0 2.6])
                        yTickValues = 0:1:2;
                    elseif iType == 1 && iRatio == 2
                        ylim([0 0.54])
                        yTickValues = 0:0.1:0.5;
                    elseif iType == 2 && iRatio == 2
                        ylim([0 0.54])
                        yTickValues = 0:0.1:0.5;
                    end
                % Upper mesopelagic
                elseif (iSubplot > 6 && iSubplot <= 12)
                    if iType == 1 && iRatio == 1
                        ylim([0 1.6])
                        yTickValues = 0:0.5:1.5;
                    elseif iType == 2 && iRatio == 1
                        ylim([0 2.6])
                        yTickValues = 0:1:2;
                    elseif iType == 1 && iRatio == 2
                        ylim([0 0.54])
                        yTickValues = 0:0.1:0.5;
                    elseif iType == 2 && iRatio == 2
                        ylim([0 0.54])
                        yTickValues = 0:0.1:0.5;
                    end
                % Lower mesopelagic
                elseif (iSubplot > 12 && iSubplot <= 18)
                    if iType == 1 && iRatio == 1
                        ylim([0 1.6])
                        yTickValues = 0:0.5:1.5;
                    elseif iType == 2 && iRatio == 1
                        ylim([0 2.6])
                        yTickValues = 0:1:2;
                    elseif iType == 1 && iRatio == 2
                        ylim([0 0.54])
                        yTickValues = 0:0.1:0.5;
                    elseif iType == 2 && iRatio == 2
                        ylim([0 0.54])
                        yTickValues = 0:0.1:0.5;
                    end
                % Base of mesopelagic
                else
                    if iType == 1 && iRatio == 1
                        ylim([0 7])
                        yTickValues = 0:2:6;
                    elseif iType == 2 && iRatio == 1
                        ylim([0 7])
                        yTickValues = 0:2:6;
                    elseif iType == 1 && iRatio == 2
                        ylim([0 1.3])
                        yTickValues = 0:0.5:1;
                    elseif iType == 2 && iRatio == 2
                        ylim([0 1.3])
                        yTickValues = 0:0.5:1;
                    end
                end
                
                yticks(yTickValues)
                yticklabels(yTickValues)

                yl = ylim; % store [ymin ymax]
                if yl(2) <= 1.9
                    ytickformat('%.1f')
                    yt = yticks;
                    ylab = arrayfun(@(x) sprintf('%.1f',x), yt, 'UniformOutput', false);
                    ylab(strcmp(ylab,'0.0')) = {'0'};
                    yticklabels(ylab);
                else
                    ytickformat('%.0f')
                end
    
                if (iSubplot == 1 || iSubplot == 7 || iSubplot == 13 || iSubplot == 19)
                    if iRatio == 1
                        yl = ylabel('PIC:POC (mol ratio)');
                    elseif iRatio == 2
                        yl = ylabel('Si:POC (mol ratio)');
                    end
                    yl.Position(1) = yl.Position(1) - 0.5;
                end
    
                set(gca, 'FontSize', 12);
    
                xlim([0.5 12+0.5])
                xticks(1:12);
                xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
                xtickangle(0); % Keep x-tick labels horizontal
                axh = gca;
                axh.XAxis.FontSize = 10; 
                
                if (iSubplot >= 1 && iSubplot <= 6)
                    tl = title(desiredLocationOrder(iLoc),'FontSize',14);
                    tl.Visible = 'on';
                    %tl.Position(2) = tl.Position(2) + 0.20;
                end
        
                grid on;
                axh.XGrid = 'off';
                axh.YGrid = 'on';
        
            end % iLoc 
        end % iDh
        
        if iRatio == 1
            ratioTag = 'pictopoc';
        elseif iRatio == 2
            ratioTag = 'bsitopoc';
        end
        if iType == 1
            typeTag = 'obs';
        elseif iType == 2
            typeTag = 'mod';
        end
        saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_',typeTag,'_',ratioTag))

    end % iType
end % iRatio

end % plotMonthlyMineralFluxRatios

% *************************************************************************

function plotHeatmapOfLocalPnumBySizeClassAndMonth(isModelled,figureSubfolderName,...
    figurePrefixName,pnumVar,pnumDepths,particleSizeClassEdges,filenameTimeseriesInformation)

nSizeClasses = length(particleSizeClassEdges)-1;

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES','STATION_TAGS')
nLocs = length(STATION_NAMES);

% Mapping properties
yaxisMax = [2000, 2000, 2000, 2000, 2000, 1000];
myColourmap = brewermap(1000,'*Spectral'); % Spectral is the equivalent of Jet

% Modify the tick labels for the colour bar
c = [0.001, 0.01, 0.1, 1, 10, 100, 1000, 1e4, 1e5, 1e6, 1e7];
cbLabels = cell(size(c)); 
for i = 1:length(c)
    if c(i) > 100
        exponent = log10(c(i));
        cbLabels{i} = ['10^' num2str(exponent)];
    else
        cbLabels{i} = num2str(c(i));
    end
end

for iLoc = 1:nLocs
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.50],'Color','w')
    haxis = zeros(nSizeClasses,1);

    for iSc = 1:nSizeClasses
        
        haxis(iSc) = subaxis(4,4,iSc,'Spacing',0.05,'Padding',0.006,'Margin',0.08);
        ax(iSc).pos = get(haxis(iSc),'Position');
        if (iSc == 1 || iSc == 5 || iSc == 9 || iSc == 13)
            ax(iSc).pos(1) = ax(iSc).pos(1);
        elseif (iSc == 2 || iSc == 6 || iSc == 10 || iSc == 14) 
            ax(iSc).pos(1) = ax(iSc).pos(1) - 0.025;
        elseif (iSc == 3 || iSc == 7 || iSc == 11 || iSc == 15)
            ax(iSc).pos(1) = ax(iSc).pos(1) - 0.050;
        elseif (iSc == 4 || iSc == 8 || iSc == 12 || iSc == 16)
            ax(iSc).pos(1) = ax(iSc).pos(1) - 0.075;   
        end
        ax(iSc).pos(2) = ax(iSc).pos(2) - 0.02;
        set(haxis(iSc),'Position',ax(iSc).pos) 

        % Extract data
        thisPnum = squeeze(pnumVar(iSc,:,:,iLoc)); 
        thisPnum(thisPnum==0) = NaN;
        theseDepths = pnumDepths(:,iLoc);

        % Discard last sampled depth in modelled data (sea floor)
        if isModelled  
            for iMonth = 1:12
                monthlyData = thisPnum(:,iMonth);
                lastIdx = find(~isnan(monthlyData),1,'last');
                if ~isempty(lastIdx)
                    thisPnum(lastIdx,iMonth) = NaN;
                end
            end
        end
    
        % Identify valid depth indices (exclude zeros)
        validDepthIdx = theseDepths ~= 0;
        y = theseDepths(validDepthIdx)';
        x = 1:12;
        dataMatrix = log10(thisPnum(validDepthIdx,:));

        % Plot
        h = imagesc(x,y,dataMatrix);
        set(h, 'AlphaData', ~isnan(dataMatrix));  % makes NaNs fully transparent
        colormap(myColourmap)
        caxis([log10(c(1)), log10(c(end))]);
        shading flat; % smooth the color transitions
        hold on
        
        % Manually add vertical grid lines
        for iCol = 1.5 : 1 : 11.5
            xline(iCol, 'w-', 'LineWidth', 1);
        end
        hold off;

        % xticks
        xticks(1:12)
        xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
        axh = gca;
        axh.XAxis.FontSize = 8;
        set(gca, 'XTickLabelRotation', 0);  % Keep x-tick labels horizontal

        % yticks
        ylim([0 yaxisMax(iLoc)])
        numTicks = 5;
        L = get(gca,'ylim');
        set(gca,'ytick',linspace(L(1),L(2),numTicks)) 

        if (iSc == 1 || iSc == 5 || iSc == 9 || iSc == 13)
            yticksValues = get(gca,'ytick');
            yticklabels(arrayfun(@num2str, yticksValues, 'UniformOutput', false));
            ytickformat(gca, '%.0f'); 
        else
            yticklabels([]);
        end
  
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.02, 0.02]; % make tick marks longer

        set(gca,'YDir','Reverse','XAxisLocation','Bottom')

        % Box: define the rectangle position and size
        xLimits = xlim;
        yLimits = ylim;
        rectangle('Position', [xLimits(1), yLimits(1), diff(xLimits), diff(yLimits)], ...
            'EdgeColor', 'k', 'LineWidth', 1);

        % Title
        if (iSc < nSizeClasses)
            title([num2str(particleSizeClassEdges(iSc)) '—' num2str(particleSizeClassEdges(iSc+1)) ' \mum']);
        else
            title(['>' num2str(particleSizeClassEdges(iSc)) ' \mum']);
        end
        
        if (iSc == 1 || iSc == 5 || iSc == 9 || iSc == 13)
            ylh = ylabel('Depth (m)');
            ylh.Position(1) = ylh.Position(1) - 0.50; 
            % ylh.Position(2) = ylh.Position(2) + 2700; 
        end
        
        if (iSc == 4)
            cb = colorbar('Location','eastoutside');
            cb.YTick = log10(c);
            cb.YTickLabel = cbLabels;
            cb.Position(1) = cb.Position(1) + 0.07;
            cb.Position(2) = cb.Position(2) - 0.67;
            cb.Position(3) = 0.015; % WIDTH
            cb.Position(4) = 0.845; % LENGTH
            cb.Label.String = 'Particle concentration (#/L)';
            cb.FontSize = 11;
        end

        set(gca,'ColorScale','log');

    end % iSc
    
    % Give common title to the figure
    a = axes; % create a new axis
    t = title(STATION_NAMES(iLoc),'FontSize',16);
    a.Visible = 'off';
    t.Visible = 'on';
    t.Position(1) = t.Position(1) - 0.07; t.Position(2) = t.Position(2) + 0.035;

    saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_',STATION_TAGS{iLoc}))
    
end % iLoc

end % plotHeatmapOfLocalPnumBySizeClassAndMonth

% *************************************************************************

function plotMonthlyRelativeContributionEsdToPocFlux(config,filenameTimeseriesInformation,...
    figureSubfolderName,figurePrefixName,fluxVar,pnumVar)

% Size class definitions
particleSizeClassEdges = config.particleDiameterClassBounds;
nSizeClasses = length(particleSizeClassEdges)-1;

% Depth definitions
fluxDepths = config.availSedTrapDeployDepths;
pnumDepths = config.availImagSysDeployDepths;

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES','STATION_TAGS')
nLocs = length(STATION_NAMES);

% Calculate contribution of each particle size class to total particle
% number by month and location
nImagSysDepths = length(pnumDepths);
pnumRelativeContribution = NaN(nSizeClasses,nImagSysDepths,12,nLocs);
for iLoc = 1:nLocs
    for iMonth = 1:12
        for iDepth = 1:nImagSysDepths
            for iSc = 1:nSizeClasses
                pnumRelativeContribution(iSc,iDepth,iMonth,iLoc) =...
                    pnumVar(iSc,iDepth,iMonth,iLoc)/sum(pnumVar(:,iDepth,iMonth,iLoc),1);
            end
        end
    end
end

% Colour choice
coloursSc = hsv(nSizeClasses);

% Legend labels
legendLabels = cell(1,nSizeClasses);
for iSc = 1:nSizeClasses
    if iSc < 16
        legendLabels{iSc} = [num2str(particleSizeClassEdges(iSc)) '–' num2str(particleSizeClassEdges(iSc+1)) ' \mum'];
    else
        legendLabels{iSc} = ['>' num2str(particleSizeClassEdges(iSc)) ' \mum'];
    end
end

% Desired depths to reveal
desiredDepths = [50,100,200,500,1000,1500,2000];
maxDesiredDepth = desiredDepths(end);

monthLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for iLoc = 1:nLocs

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.75],'Color','w')
    haxis = zeros(12,1);
    
    for iMonth = 1:12
        haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.055,'Padding',0,'Margin',0.11);
        ax(iMonth).pos = get(haxis(iMonth),'Position');
        ax(iMonth).pos(2) = ax(iMonth).pos(2) + 0.050; % move up
        ax(iMonth).pos(1) = ax(iMonth).pos(1) + 0.020; % move to the right
        set(haxis(iMonth),'Position',ax(iMonth).pos) 
    
        % .................................................................
    
        % Plot relative contributions of the pnum by size classes

        % Extract data
        thisPnum = squeeze(pnumVar(:,:,iMonth,iLoc)); 

        % Discard the last sampled depth (seafloor)
        zeroCols = all(thisPnum == 0, 1); % logical 1×48, true where all values in the column are zero
        firstZeroCol = find(zeroCols, 1); % find the first such column
        if ~isempty(firstZeroCol) && firstZeroCol > 1
            thisPnum(:,firstZeroCol-1) = NaN; 
        end

        % Limit to maximum desired depth
        [~,idxMaxDepth] = min(abs(pnumDepths - maxDesiredDepth));
        actualPnumDepths = pnumDepths(1:idxMaxDepth);
    
        % Generate normalised depth axis
        normPnumDepth = linspace(1, idxMaxDepth, idxMaxDepth);
    
        % Reshape for plotting
        [nr,nc] = size(pnumRelativeContribution(:,1:idxMaxDepth,iMonth,iLoc));
        reshRelContrib = reshape(pnumRelativeContribution(:,1:idxMaxDepth,iMonth,iLoc)',[nc,nr]);
        reshRelContrib(reshRelContrib==0) = NaN;
    
        % Plot
        hAx(1) = gca;
        hbarh = barh(hAx(1),normPnumDepth,reshRelContrib(:,:),'stacked');
        hold on
        for iSc = 1:nSizeClasses
            hbarh(iSc).FaceColor = coloursSc(iSc,:);
        end
    
        xlim(hAx(1), [0 1])
        % xticks(hAx(1), [0 0.5 1])
        % xticklabels(hAx(1), {'0','0.5','1'})
        xticklabels(hAx(1), {}) % empty xtick labels
        % xtickformat(hAx(1), '%.1f')
        hAx(1).XAxisLocation = 'top';
        hAx(1).YAxisLocation = 'left';
        hAx(1).YDir = 'reverse';
    
        % Adjust y-axis limits to remove extra space above and below the bars
        ylim(hAx(1), [min(normPnumDepth)-0.5, max(normPnumDepth)+0.5]);
    
        % Calculate y-tick positions based on closest depths
        [~,idxDesiredDepths] = arrayfun(@(d) min(abs(actualPnumDepths - d)), desiredDepths);
        ytickPositions = normPnumDepth(idxDesiredDepths);
        yticks(hAx(1), ytickPositions);

        if (iMonth == 1 || iMonth == 4 || iMonth ==7 || iMonth == 10)
            yticklabels(hAx(1), arrayfun(@num2str, desiredDepths, 'UniformOutput', false));
        else
            yticklabels([])
        end
        
        % if (iMonth == 2)
        %     xlh = xlabel('Fractional contribution of ESD class to POC flux','FontSize',12);
        %     xlh.Position(1) = xlh.Position(1) + 0.4; 
        %     xlh.Position(2) = xlh.Position(2); 
        % end
    
        % .................................................................
    
        % Plot POC flux on top
    
        hAx(1) = gca;
        hAx(2) = axes('Position',hAx(1).Position,'XAxisLocation','bottom','YAxisLocation','left',...
            'color','none','YLabel',[],'YTick',[]); % color = none, Your ax2 is using the default background color, which is white. It is hiding the other axes. 
        hold(hAx(2),'on')
    
        modelFlux = squeeze(fluxVar(:,iMonth,iLoc,1)); % mg C m-2 d-1
        modelFlux(modelFlux==0) = NaN;
    
        % Loop over each imaging depth and find the closest flux depth
        actualFluxDepths = fluxDepths;
        closestFluxDepths = zeros(size(actualPnumDepths));
        closestFlux = zeros(size(actualPnumDepths));
        for i = 1:length(actualPnumDepths)
            [~, idxDesiredDepths] = min(abs(actualFluxDepths - actualPnumDepths(i)));
            closestFluxDepths(i) = actualFluxDepths(idxDesiredDepths);
            closestFlux(i) = modelFlux(idxDesiredDepths);
        end
    
        % Generate the normalised depth axis
        normFluxDepth = linspace(1, length(closestFluxDepths), length(closestFluxDepths));
    
        % Plot
        plot(hAx(2),closestFlux,normFluxDepth,'Color','k','LineWidth',3.0,'HandleVisibility','off')
    
        % Set the y-ticks manually to ensure alignment
        ylim(hAx(2), [min(normPnumDepth)-0.5, max(normPnumDepth)+0.5]);
        yticks(hAx(1), ytickPositions);
        yticks(hAx(2), ytickPositions);
        yticklabels(hAx(2),[]); % no labels for depths of POC flux
        hAx(2).YDir = 'reverse';
    
        % Deal with x-axis labelling    
        xlim(hAx(2),[0 160])
        xTickValues = 0:75:150;
        if iLoc == 1
            xlim(hAx(2),[0 400])
            xTickValues = 0:150:300;
        elseif iLoc == 2
            xlim(hAx(2),[0 270])
            xTickValues = 0:100:200;
        elseif iLoc == 3
            xlim(hAx(2),[0 520])
            xTickValues = 0:200:400;
        elseif iLoc == 4
            xlim(hAx(2),[0 210])
            xTickValues = 0:100:200;
        end

        xticks(hAx(2),xTickValues)
        xticklabels(hAx(2),xTickValues)
        xtickformat(hAx(2),'%.0f')

        subtl = title(monthLabel(iMonth),'FontSize',16);
        subtl.Position(2) = subtl.Position(2) - 1;

    end % iMonth

    % Legend
    lg = legend(hbarh,legendLabels);
    title(lg, 'ESD classes');
    lg.Position(1) = 0.39; lg.Position(2) = -0.075;
    lg.ItemTokenSize = [11,5];
    lg.FontSize = 11;
    lg.NumColumns = 4;
    set(lg,'Box','on');

    % Give common xlabel, ylabel and title to your figure
    % Create a new axis
    a = axes;
    spt = title(STATION_NAMES(iLoc),'FontSize',18);
    xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',18);
    yl = ylabel('Depth (m)','FontSize',18);
    % Specify visibility of the current axis as 'off'
    a.Visible = 'off';
    % Specify visibility of Title, XLabel, and YLabel as 'on'
    spt.Visible = 'on';
    xl.Visible = 'on';
    yl.Visible = 'on';
    yl.Position(1) = yl.Position(1) - 0.04; yl.Position(2) = yl.Position(2) + 0.02; 
    xl.Position(1) = 0.52; xl.Position(2) = xl.Position(2) + 0.06;
    spt.Position(1) = spt.Position(1) + 0.02; spt.Position(2) = spt.Position(2) + 0.060;

    saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_',STATION_TAGS{iLoc}))

end % iLoc

end % plotMonthlyRelativeContributionEsdToPocFlux

% *************************************************************************

function plotSeasonalRelativeContributionEsdToPocFlux(config,...
    filenameTimeseriesInformation,figureSubfolderName,figurePrefixName,...
    fluxVarSeasonal,pnumVarSeasonal,pnumVarMonthly)

seasonLabel = {'Winter','Spring','Summer','Autumn'};

% Size class definitions
particleSizeClassEdges = config.particleDiameterClassBounds;
nSizeClasses = length(particleSizeClassEdges)-1;

% Depth definitions
fluxDepths = config.availSedTrapDeployDepths;
pnumDepths = config.availImagSysDeployDepths;

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')

% Rearrange model locations to match observational order
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
currentModLocationOrder = STATION_NAMES;
[~,reorderModLocIdx] = ismember(desiredLocationOrder,currentModLocationOrder); % get reordering indices
nLocs = length(STATION_NAMES);

% Calculate contribution of each particle size class to total particle
% number by month and location
nImagSysDepths = length(pnumDepths);
pnumRelativeContributionMonthly = NaN(nSizeClasses,nImagSysDepths,12,nLocs);
for iLoc = 1:nLocs
    for iMonth = 1:12
        for iDepth = 1:nImagSysDepths
            for iSc = 1:nSizeClasses
                pnumRelativeContributionMonthly(iSc,iDepth,iMonth,iLoc) =...
                    pnumVarMonthly(iSc,iDepth,iMonth,iLoc)/sum(pnumVarMonthly(:,iDepth,iMonth,iLoc),1);
            end
        end
    end
end

% Calculate the seasonal average
pnumRelativeContributionSeasonal = NaN(nSizeClasses,nImagSysDepths,4,nLocs);
for iSeason = 1:4
    pnumRelativeContributionSeasonal(:,:,iSeason,:) = mean(pnumRelativeContributionMonthly(:,:,config.seasonIndices{iSeason},:),3,'omitnan');
end

% Colour choice
coloursSc = flipud(brewermap(nSizeClasses,'Spectral'));

% Legend labels
legendLabels = cell(1,nSizeClasses);
for iSc = 1:nSizeClasses
    if iSc < 16
        legendLabels{iSc} = [num2str(particleSizeClassEdges(iSc)) '–' num2str(particleSizeClassEdges(iSc+1))];
    else
        legendLabels{iSc} = ['>' num2str(particleSizeClassEdges(iSc))];
    end
end

% Desired depths to reveal
desiredDepths = [0,100,500,1000,1500];
maxDesiredDepth = desiredDepths(end);

% Apply the right index location order
fluxVarSeasonal = fluxVarSeasonal(:,:,reorderModLocIdx,:);
pnumVarSeasonal = pnumVarSeasonal(:,:,:,reorderModLocIdx);
pnumRelativeContributionSeasonal = pnumRelativeContributionSeasonal(:,:,:,reorderModLocIdx);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.50],'Color','w')
haxis = zeros(nLocs*4,1);

for iSeason = 1:4
    for iLoc = 1:nLocs

        % Compute linear subplot index (row-wise fill)
        iSubplot = (iSeason - 1) * 6 + iLoc;
        haxis(iSubplot) = subaxis(4,nLocs,iSubplot,'Spacing',0.038,'Padding',0,'Margin',0.11);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        % Shift all plots a little bit up and to the left
        if (iSubplot <= 6)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.040;   
        elseif (iSubplot > 6 && iSubplot <= 12)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.015;        
        elseif (iSubplot > 12 && iSubplot <= 18)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.010; 
        elseif (iSubplot > 18)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.035;
        end

        if (iSubplot == 1 || iSubplot == 7 || iSubplot == 13 || iSubplot == 19)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.035;
        elseif (iSubplot == 2 || iSubplot == 8 || iSubplot == 14 || iSubplot == 20)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.041;
        elseif (iSubplot == 3 || iSubplot == 9 || iSubplot == 15 || iSubplot == 21)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.047;
        elseif (iSubplot == 4 || iSubplot == 10 || iSubplot == 16 || iSubplot == 22)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.053;
        elseif (iSubplot == 5 || iSubplot == 11 || iSubplot == 17 || iSubplot == 23)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.059;
        elseif (iSubplot == 6 || iSubplot == 12 || iSubplot == 18 || iSubplot == 24)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.065;
        end

        set(haxis(iSubplot),'Position',ax(iSubplot).pos)
    
        % .................................................................
    
        % Plot relative contributions of the pnum by size classes

        % Extract data
        thisPnum = squeeze(pnumVarSeasonal(:,:,iSeason,iLoc)); 

        % Discard the last sampled depth (seafloor)
        zeroCols = all(thisPnum == 0, 1); % logical 1×48, true where all values in the column are zero
        firstZeroCol = find(zeroCols, 1); % find the first such column
        if ~isempty(firstZeroCol) && firstZeroCol > 1
            thisPnum(:,firstZeroCol-1) = NaN; 
        end

        % Limit to maximum desired depth
        [~,idxMaxDepth] = min(abs(pnumDepths - maxDesiredDepth));
        actualPnumDepths = pnumDepths(1:idxMaxDepth);
    
        % Generate normalised depth axis
        normPnumDepth = linspace(1, idxMaxDepth, idxMaxDepth);
    
        % Reshape for plotting
        [nr,nc] = size(pnumRelativeContributionSeasonal(:,1:idxMaxDepth,iSeason,iLoc));
        reshRelContrib = reshape(pnumRelativeContributionSeasonal(:,1:idxMaxDepth,iSeason,iLoc)',[nc,nr]);
        reshRelContrib(reshRelContrib==0) = NaN;
    
        % Plot
        hAx(1) = gca;
        hbarh = barh(hAx(1),normPnumDepth,reshRelContrib(:,:),'stacked');
        hold on
        for iSc = 1:nSizeClasses
            hbarh(iSc).FaceColor = coloursSc(iSc,:);
        end
        set(hbarh, 'EdgeColor', 'k', 'LineWidth', 0.25); % makes bar plot lines rhinner
    
        xlim(hAx(1), [0 1])
        % xticks(hAx(1), [0 0.5 1])
        % xticklabels(hAx(1), {'0','0.5','1'})
        xticklabels(hAx(1), {}) % empty xtick labels
        % xtickformat(hAx(1), '%.1f')
        hAx(1).XAxisLocation = 'top';
        hAx(1).YAxisLocation = 'left';
        hAx(1).YDir = 'reverse';
    
        % Adjust y-axis limits to remove extra space above and below the bars
        ylim(hAx(1), [min(normPnumDepth)-0.5, max(normPnumDepth)+0.5]);
    
        % Calculate y-tick positions based on closest depths
        [~,idxDesiredDepths] = arrayfun(@(d) min(abs(actualPnumDepths - d)), desiredDepths);
        ytickPositions = normPnumDepth(idxDesiredDepths);
        yticks(hAx(1), ytickPositions);

        hAx(1).YAxis.TickDirection = 'out';
        hAx(1).TickLength = [0.03, 0.03]; % make tick marks longer
    
        if (iLoc == 1)
            yticklabels(hAx(1), arrayfun(@num2str, desiredDepths, 'UniformOutput', false));
        else
            yticklabels([])
        end

        % .................................................................
    
        % Plot POC flux on top
    
        hAx(1) = gca;
        hAx(2) = axes('Position',hAx(1).Position,'XAxisLocation','bottom','YAxisLocation','left',...
            'color','none','YLabel',[],'YTick',[]); % color = none, Your ax2 is using the default background color, which is white. It is hiding the other axes. 
        hold(hAx(2),'on')
    
        modelFlux = squeeze(fluxVarSeasonal(:,iSeason,iLoc,1)); % mg C m-2 d-1
        modelFlux(modelFlux==0) = NaN;
    
        % Loop over each imaging depth and find the closest flux depth
        actualFluxDepths = fluxDepths;
        closestFluxDepths = zeros(size(actualPnumDepths));
        closestFlux = zeros(size(actualPnumDepths));
        for i = 1:length(actualPnumDepths)
            [~, idxDesiredDepths] = min(abs(actualFluxDepths - actualPnumDepths(i)));
            closestFluxDepths(i) = actualFluxDepths(idxDesiredDepths);
            closestFlux(i) = modelFlux(idxDesiredDepths);
        end
    
        % Generate the normalised depth axis
        normFluxDepth = linspace(1, length(closestFluxDepths), length(closestFluxDepths));
    
        % Plot
        plot(hAx(2),closestFlux,normFluxDepth,'Color','k','LineWidth',3.0,'HandleVisibility','off')
    
        % Set the y-ticks manually to ensure alignment
        ylim(hAx(2), [min(normPnumDepth)-0.5, max(normPnumDepth)+0.5]);
        yticks(hAx(1), ytickPositions);
        yticks(hAx(2), ytickPositions);
        yticklabels(hAx(2),[]); % no labels for depths of POC flux
        hAx(2).YDir = 'reverse';
    
        % Deal with x-axis labelling    
        xlim(hAx(2),[0 120])
        xTickValues = 0:50:100;
        if iLoc == 1 % HOT/ALOHA
            xlim(hAx(2),[0 105])
            xTickValues = 0:50:100;
        elseif iLoc == 2 % BATS/OFP 
            xlim(hAx(2),[0 175])
            xTickValues = 0:75:150;
        elseif iLoc == 3 % EqPac
            xlim(hAx(2),[0 380])
            xTickValues = 0:150:300;
        elseif iLoc == 4 % PAP-SO
            xlim(hAx(2),[0 440])
            xTickValues = 0:200:400;       
        elseif iLoc == 5 % OSP
            xlim(hAx(2),[0 280])
            xTickValues = 0:100:200;
        elseif iLoc == 6
            xlim(hAx(2),[0 110])
            xTickValues = 0:50:100;
        end

        xticks(hAx(2),xTickValues)
        xticklabels(hAx(2),xTickValues)
        xtickformat(hAx(2),'%.0f')
        set(hAx(2), 'XTickLabelRotation', 0);

        if iSeason == 1
            text(0.5, 1.35, desiredLocationOrder{iLoc}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 14, 'FontWeight', 'bold');
        end

    end % iLoc
end % iSeason

% Add season labels above each row
for iSeason = 1:4
    % Get subplot indices for the current row
    rowIdx = (iSeason - 1)*6 + (1:6);

    % Extract the axis positions for the row
    rowPositions = arrayfun(@(h) get(haxis(h), 'Position'), rowIdx, 'UniformOutput', false);
    rowPositions = cat(1, rowPositions{:});  % Convert to numeric matrix: 6x4 [x y w h]

    % Compute average x center of the row (between subplots 3 and 4 is safe)
    xCenters = rowPositions(:,1) + rowPositions(:,3)/2;
    xCenterRow = mean(xCenters);

    % Compute highest y position in that row (top of the tallest axis)
    topEdges = rowPositions(:,2) + rowPositions(:,4);
    yTopRow = max(topEdges) + 0.008;  % Add small offset above the top

    % Add season label annotation
    annotation('textbox', [xCenterRow - 0.05, yTopRow, 0.1, 0.03], ...
        'String', seasonLabel{iSeason}, ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center');
end

% Legend
lg = legend(hbarh, legendLabels);
lg.Title.String = 'ESD classes (\mum)';
lg.Position(1) = 0.84;
lg.Position(2) = 0.535;
lg.ItemTokenSize = [11,4];
lg.FontSize = 11;
lg.NumColumns = 1;
set(lg, 'Box', 'on');

% Give common xlabel, ylabel and title to your figure
% Create a new axis
a = axes;
xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',14);
yl = ylabel('Depth (m)','FontSize',14);
% Specify visibility of the current axis as 'off'
a.Visible = 'off';
% Specify visibility of Title, XLabel, and YLabel as 'on'
xl.Visible = 'on';
yl.Visible = 'on';
yl.Position(1) = yl.Position(1) - 0.10; 
xl.Position(1) = 0.42; xl.Position(2) = xl.Position(2) - 0.050;

% Save
saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName))

end % plotSeasonalRelativeContributionEsdToPocFlux

% *************************************************************************

function plotSeasonalRelativeContributionVsinkToPocFlux(config,...
    filenameTimeseriesInformation,figureSubfolderName,figurePrefixName,...
    fluxVarSeasonal,pnumVarSeasonal,pnumVarMonthly)
%%
seasonLabel = {'Winter','Spring','Summer','Autumn'};

% Velocity class definitions
particleVelocityClassEdges = config.particleVelocityClassBounds;
nVeloClasses = length(particleVelocityClassEdges)-1;

% Depth definitions
fluxDepths = config.availSedTrapDeployDepths;
pnumDepths = config.availImagSysDeployDepths;

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')

% Rearrange model locations to match observational order
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
currentModLocationOrder = STATION_NAMES;
[~,reorderModLocIdx] = ismember(desiredLocationOrder,currentModLocationOrder); % get reordering indices
nLocs = length(STATION_NAMES);

% Calculate contribution of each particle velocity class to total particle
% number by month and location
nImagSysDepths = length(pnumDepths);
pnumRelativeContributionMonthly = NaN(nVeloClasses,nImagSysDepths,12,nLocs);
for iLoc = 1:nLocs
    for iMonth = 1:12
        for iDepth = 1:nImagSysDepths
            for iVc = 1:nVeloClasses
                pnumRelativeContributionMonthly(iVc,iDepth,iMonth,iLoc) =...
                    pnumVarMonthly(iVc,iDepth,iMonth,iLoc)/sum(pnumVarMonthly(:,iDepth,iMonth,iLoc),1);
            end
        end
    end
end

% Calculate the seasonal average
pnumRelativeContributionSeasonal = NaN(nVeloClasses,nImagSysDepths,4,nLocs);
for iSeason = 1:4
    pnumRelativeContributionSeasonal(:,:,iSeason,:) = mean(pnumRelativeContributionMonthly(:,:,config.seasonIndices{iSeason},:),3,'omitnan');
end

% Colour choice
coloursVc = flipud(brewermap(nVeloClasses,'Spectral'));

% Legend labels
legendLabels = cell(1,nVeloClasses);
for iVc = 1:nVeloClasses
    if iVc < 16
        % Format lower and upper edges
        lowEdge  = particleVelocityClassEdges(iVc);
        highEdge = particleVelocityClassEdges(iVc+1);

        if lowEdge < 1
            lowStr = sprintf('%.1f', lowEdge);   % two decimals
        else
            lowStr = sprintf('%.0f', lowEdge);   % no decimals
        end

        if highEdge < 1
            highStr = sprintf('%.1f', highEdge);
        else
            highStr = sprintf('%.0f', highEdge);
        end

        legendLabels{iVc} = [lowStr '–' highStr];
    else
        edgeVal = particleVelocityClassEdges(iVc);
        if edgeVal < 1
            edgeStr = sprintf('%.1f', edgeVal);
        else
            edgeStr = sprintf('%.0f', edgeVal);
        end
        legendLabels{iVc} = ['>' edgeStr];
    end
end

% Desired depths to reveal
desiredDepths = [0,100,500,1000,1500];
maxDesiredDepth = desiredDepths(end);

% Apply the right index location order
fluxVarSeasonal = fluxVarSeasonal(:,:,reorderModLocIdx,:);
pnumVarSeasonal = pnumVarSeasonal(:,:,:,reorderModLocIdx);
pnumRelativeContributionSeasonal = pnumRelativeContributionSeasonal(:,:,:,reorderModLocIdx);
%%
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.50],'Color','w')
haxis = zeros(nLocs*4,1);

for iSeason = 1:4
    for iLoc = 1:nLocs

        % Compute linear subplot index (row-wise fill)
        iSubplot = (iSeason - 1) * 6 + iLoc;
        haxis(iSubplot) = subaxis(4,nLocs,iSubplot,'Spacing',0.038,'Padding',0,'Margin',0.11);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        % Shift all plots a little bit up and to the left
        if (iSubplot <= 6)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.040;   
        elseif (iSubplot > 6 && iSubplot <= 12)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.015;        
        elseif (iSubplot > 12 && iSubplot <= 18)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.010; 
        elseif (iSubplot > 18)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.035;
        end

        if (iSubplot == 1 || iSubplot == 7 || iSubplot == 13 || iSubplot == 19)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.035;
        elseif (iSubplot == 2 || iSubplot == 8 || iSubplot == 14 || iSubplot == 20)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.041;
        elseif (iSubplot == 3 || iSubplot == 9 || iSubplot == 15 || iSubplot == 21)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.047;
        elseif (iSubplot == 4 || iSubplot == 10 || iSubplot == 16 || iSubplot == 22)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.053;
        elseif (iSubplot == 5 || iSubplot == 11 || iSubplot == 17 || iSubplot == 23)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.059;
        elseif (iSubplot == 6 || iSubplot == 12 || iSubplot == 18 || iSubplot == 24)
            ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.065;
        end

        set(haxis(iSubplot),'Position',ax(iSubplot).pos)
    
        % .................................................................
    
        % Plot relative contributions of the pnum by size classes

        % Extract data
        thisPnum = squeeze(pnumVarSeasonal(:,:,iSeason,iLoc)); 

        % Discard the last sampled depth (seafloor)
        zeroCols = all(thisPnum == 0, 1); % logical 1×48, true where all values in the column are zero
        firstZeroCol = find(zeroCols, 1); % find the first such column
        if ~isempty(firstZeroCol) && firstZeroCol > 1
            thisPnum(:,firstZeroCol-1) = NaN; 
        end

        % Limit to maximum desired depth
        [~,idxMaxDepth] = min(abs(pnumDepths - maxDesiredDepth));
        actualPnumDepths = pnumDepths(1:idxMaxDepth);
    
        % Generate normalised depth axis
        normPnumDepth = linspace(1, idxMaxDepth, idxMaxDepth);
    
        % Reshape for plotting
        [nr,nc] = size(pnumRelativeContributionSeasonal(:,1:idxMaxDepth,iSeason,iLoc));
        reshRelContrib = reshape(pnumRelativeContributionSeasonal(:,1:idxMaxDepth,iSeason,iLoc)',[nc,nr]);
        reshRelContrib(reshRelContrib==0) = NaN;

        % Plot
        hAx(1) = gca;
        hbarh = barh(hAx(1),normPnumDepth,reshRelContrib(:,:),'stacked');
        hold on
        for iVc = 1:nVeloClasses
            hbarh(iVc).FaceColor = coloursVc(iVc,:);
        end
        set(hbarh, 'EdgeColor', 'k', 'LineWidth', 0.25); % makes bar plot lines rhinner

        xlim(hAx(1), [0 1])
        xticklabels(hAx(1), {}) % empty xtick labels
        hAx(1).XAxisLocation = 'top';
        hAx(1).YAxisLocation = 'left';
        hAx(1).YDir = 'reverse';
    
        % Adjust y-axis limits to remove extra space above and below the bars
        ylim(hAx(1), [min(normPnumDepth)-0.5, max(normPnumDepth)+0.5]);
    
        % Calculate y-tick positions based on closest depths
        [~,idxDesiredDepths] = arrayfun(@(d) min(abs(actualPnumDepths - d)), desiredDepths);
        ytickPositions = normPnumDepth(idxDesiredDepths);
        yticks(hAx(1), ytickPositions);

        hAx(1).YAxis.TickDirection = 'out';
        hAx(1).TickLength = [0.03, 0.03]; % make tick marks longer
    
        if iLoc == 1
            yticklabels(hAx(1), arrayfun(@num2str, desiredDepths, 'UniformOutput', false));
        else
            yticklabels([])
        end

        % .................................................................
    
        % Plot POC flux on top
    
        hAx(1) = gca;
        hAx(2) = axes('Position',hAx(1).Position,'XAxisLocation','bottom','YAxisLocation','left',...
            'color','none','YLabel',[],'YTick',[]); % color = none, Your ax2 is using the default background color, which is white. It is hiding the other axes. 
        hold(hAx(2),'on')
    
        modelFlux = squeeze(fluxVarSeasonal(:,iSeason,iLoc,1)); % mg C m-2 d-1
        modelFlux(modelFlux==0) = NaN;
    
        % Loop over each imaging depth and find the closest flux depth
        actualFluxDepths = fluxDepths;
        closestFluxDepths = zeros(size(actualPnumDepths));
        closestFlux = zeros(size(actualPnumDepths));
        for i = 1:length(actualPnumDepths)
            [~, idxDesiredDepths] = min(abs(actualFluxDepths - actualPnumDepths(i)));
            closestFluxDepths(i) = actualFluxDepths(idxDesiredDepths);
            closestFlux(i) = modelFlux(idxDesiredDepths);
        end
    
        % Generate the normalised depth axis
        normFluxDepth = linspace(1, length(closestFluxDepths), length(closestFluxDepths));
    
        % Plot
        plot(hAx(2),closestFlux,normFluxDepth,'Color','k','LineWidth',3.0,'HandleVisibility','off')
    
        % Set the y-ticks manually to ensure alignment
        ylim(hAx(2), [min(normPnumDepth)-0.5, max(normPnumDepth)+0.5]);
        yticks(hAx(1), ytickPositions);
        yticks(hAx(2), ytickPositions);
        yticklabels(hAx(2),[]); % no labels for depths of POC flux
        hAx(2).YDir = 'reverse';
    
        % Deal with x-axis labelling    
        xlim(hAx(2),[0 120])
        xTickValues = 0:50:100;
        if iLoc == 1 % HOT/ALOHA
            xlim(hAx(2),[0 105])
            xTickValues = 0:50:100;
        elseif iLoc == 2 % BATS/OFP 
            xlim(hAx(2),[0 175])
            xTickValues = 0:75:150;
        elseif iLoc == 3 % EqPac
            xlim(hAx(2),[0 380])
            xTickValues = 0:150:300;
        elseif iLoc == 4 % PAP-SO
            xlim(hAx(2),[0 440])
            xTickValues = 0:200:400;       
        elseif iLoc == 5 % OSP
            xlim(hAx(2),[0 280])
            xTickValues = 0:100:200;
        elseif iLoc == 6
            xlim(hAx(2),[0 110])
            xTickValues = 0:50:100;
        end

        xticks(hAx(2),xTickValues)
        xticklabels(hAx(2),xTickValues)
        xtickformat(hAx(2),'%.0f')
        set(hAx(2), 'XTickLabelRotation', 0);

        if iSeason == 1
            text(0.5, 1.35, desiredLocationOrder{iLoc}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 14, 'FontWeight', 'bold');
        end

    end % iLoc
end % iSeason

% Add season labels above each row
for iSeason = 1:4
    % Get subplot indices for the current row
    rowIdx = (iSeason - 1)*6 + (1:6);

    % Extract the axis positions for the row
    rowPositions = arrayfun(@(h) get(haxis(h), 'Position'), rowIdx, 'UniformOutput', false);
    rowPositions = cat(1, rowPositions{:});  % Convert to numeric matrix: 6x4 [x y w h]

    % Compute average x center of the row (between subplots 3 and 4 is safe)
    xCenters = rowPositions(:,1) + rowPositions(:,3)/2;
    xCenterRow = mean(xCenters);

    % Compute highest y position in that row (top of the tallest axis)
    topEdges = rowPositions(:,2) + rowPositions(:,4);
    yTopRow = max(topEdges) + 0.008;  % Add small offset above the top

    % Add season label annotation
    annotation('textbox', [xCenterRow - 0.05, yTopRow, 0.1, 0.03], ...
        'String', seasonLabel{iSeason}, ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center');
end

% Legend
lg = legend(hbarh, legendLabels);
lg.Title.String = 'Vp classes (m/d)';
lg.Position(1) = 0.86;
lg.Position(2) = 0.535;
lg.ItemTokenSize = [11,4];
lg.FontSize = 11;
lg.NumColumns = 1;
set(lg, 'Box', 'on');

% Give common xlabel, ylabel and title to your figure
% Create a new axis
a = axes;
xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',14);
yl = ylabel('Depth (m)','FontSize',14);
% Specify visibility of the current axis as 'off'
a.Visible = 'off';
% Specify visibility of Title, XLabel, and YLabel as 'on'
xl.Visible = 'on';
yl.Visible = 'on';
yl.Position(1) = yl.Position(1) - 0.10; 
xl.Position(1) = 0.42; xl.Position(2) = xl.Position(2) - 0.050;

% Save
saveFigureInFolder(figureSubfolderName,figurePrefixName)

end % plotSeasonalRelativeContributionVsinkToPocFlux

% *************************************************************************

function plotAverageParticleAttributes(particleAttributesByScAnnualDhMetadata,...
    particleAttributesByScSeasonalDhMetadata,fluxMonthlyDh,fluxSeasonalDh,...
    figureSubfolderName,figurePrefixName)
%%
locationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
sizeLabels = {'Small (<150 \mum ESD)','Large (>150 \mum ESD)'};
depthLabels = {'Base of euphotic','Upper mesopelagic','Lower mesopelagic','Bathypelagic/mesopelagic'};

nDepthLevels = numel(depthLabels);
nAttr = length(particleAttributesByScAnnualDhMetadata);
nSizeCat = length(sizeLabels);
nLocs = length(locationOrder);

%% First figure: bar plots

colourScheme = flipud(brewermap(nLocs,'*spectral'));

for iAtt = 1:nAttr
    attributeArray = particleAttributesByScAnnualDhMetadata{iAtt,3};
    attributeNameTag = particleAttributesByScAnnualDhMetadata{iAtt,1};
    attributeNameLabel = particleAttributesByScAnnualDhMetadata{iAtt,4};

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.20 0.60],'Color','w')
    haxis = zeros(nDepthLevels,1);

    for iDh = 1:nDepthLevels
        haxis(iDh) = subaxis(nDepthLevels,1,iDh,'Spacing',0.02,'Padding',0.02,'Margin',0.06);
        ax(iDh).pos = get(haxis(iDh),'Position');
        ax(iDh).pos(1) = ax(iDh).pos(1) + 0.05;
        ax(iDh).pos(2) = ax(iDh).pos(2) + 0.05;
        set(haxis(iDh), 'Position', ax(iDh).pos)

        % Extract data and transform to logarithmic scale if necessary
        data = attributeArray(:,:,iDh)'; % 6 locations x 2 sizes --> 2 sizes x 6 locations

        % Plot
        h = bar(data); 

        % Set colours for each bar group
        for k = 1:length(h)
            h(k).FaceColor = colourScheme(k,:);
        end
        box on

        xticklabels(sizeLabels);

        if (iAtt <= 6 || iAtt == 7 || iAtt == 9 || iAtt == 10 || iAtt == 13)
            set(gca,'Yscale','log');
        end

        if iAtt == 1 % POC
            ylim([1e-14 1e-7])
            yticks([1e-12 1e-10 1e-8])
        elseif iAtt == 2 % TEP-C
            ylim([1e-18 1e-7])
        elseif iAtt == 3 % stick
            ylim([1e-5 1])
            yticks([1e-4 1e-2 1])
        elseif iAtt == 4 % CaCO3
            ylim([1e-16 1e-8])
            yticks([1e-15 1e-10])
        elseif iAtt == 5 % bSi
            ylim([1e-20 1e-7])
        elseif iAtt == 6 % clay
            ylim([1e-22 1e-8])
            yticks([1e-20 1e-10])
        elseif iAtt == 7 % volume
            ylim([1e-1 1e8])
        elseif iAtt == 8 % D3
            ylim([1 3])
        elseif iAtt == 9 % ESD
            ylim([1e0 1e4])
        elseif iAtt == 10 % rpp
            ylim([0 1e4])    
        elseif iAtt == 11 % porosity
            ylim([1e-10 1])
        elseif iAtt == 12 % rho
            ylim([0.8 2.8])
        elseif iAtt == 13 % vsink
            ylim([1e-3 10])
        elseif iAtt == 14 % pnum
            ylim([1e-2 2e8])
        end

        ylabel(attributeNameLabel);
        title(depthLabels(iDh),'FontSize',14);
        grid on

    end % iDh

    % Unique legend at the bottom
    % Create invisible axes over entire figure
    hLegendAxes = axes('Position',[0, 0, 1, 1],'Visible','off');

    % Create the legend from the bar handles of the last plot (or any subplot)
    hLegend = legend(hLegendAxes,h,locationOrder,'Orientation','horizontal','Box','off');                    
    set(hLegend, 'Units', 'normalized');
    set(hLegend, 'Position', [0.4, 0.03, 0.3, 0.05]); % [x, y, width, height]
    hLegend.NumColumns = 3;
    hLegend.ItemTokenSize = [12, 8];

    % Save
    saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_',attributeNameTag))

end % iAtt

%% Correlation matrix against biogeochemical fluxes

fluxPoc = squeeze(fluxSeasonalDh(:,:,:,1)); % 4 nz x 4 seasons x 6 locs x 3 tracers
fluxPic = squeeze(fluxSeasonalDh(:,:,:,2));
fluxBSi = squeeze(fluxSeasonalDh(:,:,:,3));
fluxMatrices = {fluxPoc,fluxPic,fluxBSi};
fluxLabels = {'POC flux','PIC flux','bSi flux'};

depthSel = 1:3; % first three depths

corrValues = nan(nAttr,nLocs,length(fluxMatrices),nSizeCat);
pValues    = nan(size(corrValues));

for iSize = 1:nSizeCat
    for iFlux = 1:length(fluxMatrices)
        for iLoc = 1:nLocs
            % Pool flux across first 3 depths
            fluxData = squeeze(fluxMatrices{iFlux}(depthSel,:,iLoc)); % 3 depths x 4 seasons
            y = fluxData(:); % from 3 x 4 --> 12
            for iAttr = 1:nAttr
                attrVal = squeeze(particleAttributesByScSeasonalDhMetadata{iAttr,3}(iLoc,iSize,depthSel,:)); % 3 depths x 4 seasons
                x = attrVal(:); % from 3 x 4 --> 12
                valid = ~isnan(x) & ~isnan(y);
                if sum(valid) > 1
                    [corrValues(iAttr,iLoc,iFlux,iSize),pValues(iAttr,iLoc,iFlux,iSize)] =...
                        corr(x(valid),y(valid),'Type','Spearman');
                else
                    [corrValues(iAttr,iLoc,iFlux,iSize),pValues(iAttr,iLoc,iFlux,iSize)] = NaN;
                end
            end
        end
    end
end

%% Second figure: correlation heatmap

sizeLabels = {'Small','Large'};

layout = struct(...
    'leftMargin', 0.11, ...
    'rightMargin', 0.14, ...
    'topMargin', 0.07, ...
    'bottomMargin', 0.22, ...
    'gapX', 0.07, ...
    'gapY', 0.03, ...
    'nCols', 1, ...
    'nRows', 1 ...
);

width = (1 - layout.leftMargin - layout.rightMargin - (layout.nCols - 1) * layout.gapX) / layout.nCols;
height = (1 - layout.topMargin - layout.bottomMargin - (layout.nRows - 1) * layout.gapY) / layout.nRows;
axHandles = gobjects(layout.nRows*layout.nCols,1);
axPositions = zeros(layout.nRows*layout.nCols,4);

% Colour bar colours
edges = -1:0.1:1;   % interval edges
nBins = length(edges)-1;  % should be 10 bins
cmap = cmocean('balance', nBins);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.30],'Color','w')

for iFlux = 1 %:length(fluxMatrices) – as in nCols
    subplotIdx = iFlux; 

    row = floor((subplotIdx-1) / layout.nCols);
    col = mod((subplotIdx-1), layout.nCols); 
    left = layout.leftMargin + col * (width + layout.gapX);
    bottom = 1 - layout.topMargin - (row + 1) * height - row * layout.gapY;
    ax = axes('Position', [left bottom width height]);
    axHandles(subplotIdx) = ax;
    axPositions(subplotIdx,:) = [left bottom width height];

    % Plot
    for iLoc = 1:nLocs
        data = squeeze(corrValues(:,:,iFlux,:)); % 10 atts x 6 locs x 2 sizes
        heatmapData = [data(:,:,1), data(:,:,2)]; % concatenate small and large side by side --> 10 x 12
        
        pval = squeeze(pValues(:,:,iFlux,:));
        heatmapPval = [pval(:,:,1), pval(:,:,2)];

        % Plot
        h = imagesc(ax,heatmapData);
    
        % Set transparency: fully opaque (1) for real numbers, transparent (0) for NaNs
        set(h, 'AlphaData', ~isnan(heatmapData));  % NaNs get 0 (transparent), others get 1
    
        colormap(ax,cmap); % colormap(ax,cmocean('balance'));
        caxis(ax,[-1 1]);

        axis tight
        % pbaspect(ax,[1.5 1 1]); % makes cells 1.5× wider than tall
    
        % Set y-ticks with latex labels
        set(ax, 'YTick', 1:nAttr, ...
            'YTickLabel', particleAttributesByScSeasonalDhMetadata(:,2), ...
            'TickLabelInterpreter', 'latex',...
            'FontSize', 12)
    
        % Set x-ticks labels
        % set(ax,'XTick',1:nLocs,'XTickLabel',locationOrder);
        xticks(1:(2*nLocs));
        xticklabels([locationOrder, locationOrder]);
        xtickangle(45);   % or 90 for vertical
    
        title(fluxLabels{iFlux},'FontSize',14);

        % Draw vertical line to mark separator
        hold on;
        xline(6.5,'k','LineWidth',2);
        hold on;

        % Add text labels (correlation values)
        hold on;
        for iAttr = 1:nAttr
            for jCol = 1:(2*nLocs)
                if ~isnan(heatmapData(iAttr,jCol))
                    
                    % Format correlation value
                    val = heatmapData(iAttr,jCol);
                    p = heatmapPval(iAttr,jCol);

                    % Choose text color: white for dark backgrounds, black otherwise
                    if abs(val) > 0.4
                        textColour = 'w';
                    else
                        textColour = 'k';
                    end

                    if p <= 0.05
                        labelStr = sprintf('*'); %sprintf('%.2f', val);
                    else
                        labelStr = '';
                    end
                    
                    % Place text at (column,row)
                    text(jCol, iAttr, labelStr, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize', 9.5, ...
                        'Color', textColour);  % you can switch 'w' for dark cells
                end
            end
        end
        hold off;

        % % Add text labels
        % threshold = 0.4;  % Set threshold for white text
        % for iRow = 1:size(heatmapData,1)
        %     for iCol = 1:size(heatmapData,2)
        %         val = data(iRow,iCol);
        %         sig = pval(iRow,iCol);
        %         % Choose text color: white for dark backgrounds, black otherwise
        %         if abs(val) > threshold
        %             textColour = 'w';
        %         else
        %             textColour = 'k';
        %         end
        % 
        %         if sig <= 0.1
        %             str = sprintf('%.2f', val);
        %         else
        %             str = '';
        %         end
        % 
        %         text(iCol, iRow, str, ...
        %             'HorizontalAlignment', 'center', ...
        %             'Color', textColour, 'FontSize', 9.5);
        %     end
        % end
    
        if iFlux == 1 %length(fluxMatrices)
            cb = colorbar(ax); % anchor to last subplot
            % cb.Ticks = edges; % ticks every 0.1
            % cb.TickLabels = arrayfun(@num2str, edges, 'UniformOutput', false); 
            
            cb.Ticks = -1:0.1:1;   % ticks every 0.1
            labels = repmat({''}, size(cb.Ticks));
            for i = 1:numel(cb.Ticks)
                if abs(round(cb.Ticks(i)/0.2)*0.2 - cb.Ticks(i)) < 1e-6
                    labels{i} = sprintf('%.1f', cb.Ticks(i));
                end
            end
            cb.TickLabels = labels;
            
            % Apply labels
            cb.TickLabels = labels;

            cb.Position(1) = cb.Position(1) + 0.09;
            cb.Label.String = 'Correlation'; 
            cb.Label.Interpreter = 'none';
            cb.FontSize = 12;               
        end

    end

end % iFlux

% Save
saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_correlation'))

end % plotAverageParticleAttributes

% *************************************************************************

function plotSinkTerms(figureSubfolderName,testrunDir,figurePrefixName,...
    config,fullpathModelRunsDir,filenameSlamsOutput,choiceTypeGridDomain)

%% SECTION 1: Load data

load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')

% Extrac annual SMS data (nIdxs x nDepths x nLocs)
annualSms = output.sms.annual;  

%% SECTION 2: Choose sink terms

% Loss terms for POC and PIC (bSi is only lost through abiotic dissolution)
idxPocLossTerms = [config.sms.idxMicrobRespOrgC,config.sms.idxMicrobRespTepC,...
    config.sms.idxZooRespOrgC,config.sms.idxZooRespTepC,...
    config.sms.idxMicrobSolubOrgC,config.sms.idxMicrobSolubTepC,...
    config.sms.idxZooExcretOrgC,config.sms.idxZooExcretTepC,...
    config.sms.idxPhotoTepC];
idxPicLossTerms = [config.sms.idxDissolCaCO3,config.sms.idxZooDissolCaCO3];

pocLossLabels = {'Microbial respiration','Mesozooplankton respiration',...
    'Microbial solubilisation','Mesozooplankton excretion','Photolysis'};
picLossLabels = {'Abiotic CaCO3 dissolution','Mesozoo. gut CaCO3 dissolution'};
nPocLossTermsTemp = length(idxPocLossTerms);
nPicLossTerms = length(idxPicLossTerms);

%% SECTION 3: Extract specific sink terms and transform into fractions

pocLossTemp = zeros(nPocLossTermsTemp,size(config.gridDepths,1),config.nLocs);
picLoss = zeros(nPicLossTerms,size(config.gridDepths,1),config.nLocs);
for iLoc = 1:config.nLocs % make sure this is the same order as output matrix
    nLocalDepths = sum(~isnan(config.gridDepths(:,iLoc)));
    for iZ = 1:nLocalDepths
        for iPoc = 1:nPocLossTermsTemp
            pocLossTemp(iPoc,iZ,iLoc) = annualSms(idxPocLossTerms(iPoc),iZ,iLoc);
        end
        for iPic = 1:nPicLossTerms
            picLoss(iPic,iZ,iLoc) = annualSms(idxPicLossTerms(iPic),iZ,iLoc);
        end
    end
end

% For POC, add orgC and TEC
nPocLossTerms = floor(nPocLossTermsTemp/2)+1;
pocLoss = zeros(nPocLossTerms,size(config.gridDepths,1),config.nLocs);
for iLoc = 1:config.nLocs % make sure this is the same order as output matrix
    nLocalDepths = sum(~isnan(config.gridDepths(:,iLoc)));
    for iZ = 1:nLocalDepths
        for iTerm = 1:nPocLossTerms
            if iTerm == 1
                iAdd = 1;
                pocLoss(iTerm,iZ,iLoc) = pocLossTemp(iAdd,iZ,iLoc) + pocLossTemp(iAdd+1,iZ,iLoc);
            elseif iTerm > 1 && iTerm < nPocLossTerms
                iAdd = iAdd + 2;
                pocLoss(iTerm,iZ,iLoc) = pocLossTemp(iAdd,iZ,iLoc) + pocLossTemp(iAdd+1,iZ,iLoc);
            elseif iTerm == nPocLossTerms
                iAdd = iAdd + 2;
                pocLoss(iTerm,iZ,iLoc) = pocLossTemp(iAdd,iZ,iLoc);
            end
        end
    end
end

% Transform quantities into fractions
fracPocLoss = zeros(nPocLossTerms,size(config.gridDepths,1),config.nLocs);
fracPicLoss = zeros(nPicLossTerms,size(config.gridDepths,1),config.nLocs);
for iLoc = 1:config.nLocs
    for iZ = 1:size(config.gridDepths,1)
        totalPocLoss = sum(pocLoss(:,iZ,iLoc));
        totalPicLoss = sum(picLoss(:,iZ,iLoc));
        for iPoc = 1:nPocLossTerms
            fracPocLoss(iPoc,iZ,iLoc) = pocLoss(iPoc,iZ,iLoc)/totalPocLoss;
        end
        for iPic = 1:nPicLossTerms
            fracPicLoss(iPic,iZ,iLoc) = picLoss(iPic,iZ,iLoc)/totalPicLoss;
        end
    end
end

%% SECTION 4: Plot

% Rearrange model locations
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'};
desiredModLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
[~,reorderModLocIdx] = ismember(desiredModLocationOrder,currentModLocationOrder); % get reordering indices

% Reorder to match desired locations
fracPocLoss = fracPocLoss(:,:,reorderModLocIdx);
fracPicLoss = fracPicLoss(:,:,reorderModLocIdx);
lossTermsDepths = config.gridDepths(:,reorderModLocIdx);

% Plot
coloursPoc = brewermap(nPocLossTerms,'Paired');
coloursPic = brewermap(nPicLossTerms,'Paired');

for iTracer = 1:2
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.65 0.25],'Color','w')
    haxis = zeros(config.nLocs,1);
    
    for iLoc = 1:config.nLocs
        haxis(iLoc) = subaxis(1,config.nLocs,iLoc,'Spacing',0.035,'Padding',0.001,'Margin',0.11);
        ax(iLoc).pos = get(haxis(iLoc),'Position');
        ax(iLoc).pos(1) = ax(iLoc).pos(1) - 0.075;
        set(haxis(iLoc),'Position',ax(iLoc).pos) 

        switch iTracer
            case 1
                y = fracPocLoss(:,:,iLoc);
                cp = coloursPoc;
                nLts = nPocLossTerms;
            case 2
                y = fracPicLoss(:,:,iLoc);
                cp = coloursPic;
                nLts = nPicLossTerms;
        end
    
        harea = area(haxis(iLoc),lossTermsDepths(:,iLoc),y'); % transpose y from nLt x nDepths --> nDepts x nLt
        for iLt = 1:nLts
            harea(iLt).FaceColor = cp(iLt,:);
        end   
        ylim([0 1])
        hold off

        set(haxis(iLoc), 'XDir', 'reverse');
        set(haxis(iLoc), 'XTick', [0 500 1000 1500 2000 3000 4000]);
        xtickformat('%d')
        if iLoc == 1
            xlabel('Depth (m)');
        end
        box on    
        title(desiredModLocationOrder(iLoc),'FontSize',14);
        view([90 -90])

    end % iLoc

    % Legend
    if iTracer == 1
        lg = legend(harea,pocLossLabels);
        lg.Position(2) = 0.61;
        lg.Title.String = 'POC sink terms';
    elseif iTracer == 2
        lg = legend(harea,picLossLabels);
        lg.Position(2) = 0.77;
        lg.Title.String = 'PIC sink terms';
    end
    lg.Position(1) = 0.82;
    lg.ItemTokenSize = [11,5];
    lg.FontSize = 10.5;
    set(lg,'Box','off');

    if (iTracer == 1)
        saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_poc'))
    elseif (iTracer == 2)
        saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_pic'))
    end

end % iTracer

end % plotSinkTerms

% *************************************************************************

function plotEncounterTerms(figureSubfolderName,testrunDir,figurePrefixName,...
    config,fullpathModelRunsDir,filenameSlamsOutput,choiceTypeGridDomain)

%% SECTION 1: Load data

load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')

% Extract annual aux. data (nIdxs x nDepths x nLocs)
annualAux = output.aux.annual;  

% Rearrange model locations
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'};
desiredModLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
[~,reorderModLocIdx] = ismember(desiredModLocationOrder,currentModLocationOrder); % get reordering indices

%% SECTION 2: Plot coagulation statistics

% Reorder to match desired locations
annualAux = annualAux(:,:,reorderModLocIdx);
nDepths = size(annualAux,2);
auxTermsDepths = config.gridDepths(1:nDepths,reorderModLocIdx);

for iFigure = 1:2
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.60 0.25],'Color','w')
    haxis = zeros(1,config.nLocs);

    for iLoc = 1:config.nLocs
        haxis(iLoc) = subaxis(1,config.nLocs,iLoc,'Spacing',0.03,'Padding',0.001,'Margin',0.08);
        ax(iLoc).pos = get(haxis(iLoc),'Position');
        ax(iLoc).pos(1) = ax(iLoc).pos(1) - 0.038;
        set(haxis(iLoc),'Position',ax(iLoc).pos)

        depths = auxTermsDepths(:,iLoc);

        if iFigure == 1
            %semilogx(haxis(iLoc),squeeze(annualAux(config.aux.idxCollisionKernel,:,iLoc))',depths,'k-','LineWidth',2); hold on;
            semilogx(haxis(iLoc),squeeze(annualAux(config.aux.idxCoagulationProbability,:,iLoc))',depths,'k-','LineWidth',2); hold on;
            semilogx(haxis(iLoc),squeeze(annualAux(config.aux.idxCoagulationSuccess,:,iLoc))',depths,'Color', [0.7 0.7 0.7],'LineWidth',2); hold off;
            xlim([1e-5 1])
            xticks([1e-5 1e-3 1e-1])
            xticklabels({'10^{-5}','10^{-3}','10^{-1}'})
            ylim([0 500])
            yticks(0:100:500)
            grid on
            set(gca,'XAxisLocation','bottom','YAxisLocation','left','Ydir','reverse','XMinorTick','off','XMinorGrid','off')
        elseif iFigure == 2
            colourPalette = brewermap(3,'Set1');
            semilogx(haxis(iLoc),squeeze(annualAux(config.aux.idxBrownianKernel,:,iLoc))',depths,'Color',colourPalette(1,:),'LineWidth',2); hold on;
            semilogx(haxis(iLoc),squeeze(annualAux(config.aux.idxShearKernel,:,iLoc))',depths,'Color',colourPalette(2,:),'LineWidth',2); hold on;
            semilogx(haxis(iLoc),squeeze(annualAux(config.aux.idxSettlingKernel,:,iLoc))',depths,'Color',colourPalette(3,:),'LineWidth',2); hold off;
            xlim([1e-16 1])
            xticks([1e-15 1e-10 1e-5 1])
            xticklabels({'10^{-15}','10^{-10}','10^{-5}','1'})
            set(gca, 'XTickLabelRotation', 0)
            ylim([0 500])
            yticks(0:100:500)
            grid on
            set(gca,'XAxisLocation','bottom','YAxisLocation','left','Ydir','reverse')
        end

        title(desiredModLocationOrder(iLoc),'FontSize',12);
        if iLoc == 1
            ylabel('Depth (m)');
        end
    
    end % iLoc

    % Legend
    if iFigure == 1
        lglabels = {'Coagu. prob.','Coagu. success'};
    elseif iFigure == 2
        lglabels = {'Brownian (m^{3} s^{-1})','Shear (m^{3} s^{-1})','Settling (m^{3} s^{-1})'};
    end
    lg = legend(haxis(6),lglabels);
    lg.Position(1) = 0.89; 
    lg.Position(2) = 0.54;
    lg.ItemTokenSize = [22,1];
    lg.FontSize = 11;
    set(lg,'Box','on');
    
    if iFigure == 1
        saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_success'))
    elseif iFigure == 2
        saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_particles'))
    end

end % iFigure

%% SECTION 3: Zooplankton encounter statistics

% Reorder to match desired locations
zooEncStats = annualAux(config.aux.idxNumParticlesEvalEncounter:config.aux.idxNumParticlesOmittedZoo,:,reorderModLocIdx);
zooEncStats(zooEncStats==0) = 0;
nDepths = size(zooEncStats,2);
auxTermsDepths = config.gridDepths(1:nDepths,reorderModLocIdx);

% Convert encounter statistics into fractional values
for iLoc = 1:config.nLocs
    locData = zooEncStats(:,:,iLoc); % 5 x 401
    reshZooEncStats = locData'; % 401 x 5
    
    % Ensure that probabilities in columns 1 and 2 sum to 1 (column 2 is
    % already a fraction of 1, whereas column 1 has total particle number).
    % Therefore, 
    reshZooEncStats(:,1) = 1 - reshZooEncStats(:,2);
    
    % Extract encounter outcome fractions (columns 3, 4 and 5 add up to 1)
    zooFrag = reshZooEncStats(:,3); 
    zooIng = reshZooEncStats(:,4);
    zooOmit = reshZooEncStats(:,5);

    % Columns 3, 4, and 5 should sum to 1. If they do not, it indicates the 
    % presence of a fourth category: particles encountered but not evaluated 
    % for ingestion. These cases occur when fragmentation evaluation fails,
    % e.g., the particle size ratio was appropriate but only a single 
    % particle was present. Such cases are reassignedto the omission pool.
    encounteredBlock = sum([zooFrag,zooOmit,zooIng],2);
    for iDh = 1:length(encounteredBlock)
        if encounteredBlock(iDh) ~= 1
            zooOmit(iDh) = zooOmit(iDh) + (1-encounteredBlock(iDh));
        end
    end

    % Consolidate adjusted fractions: [p(not-encountered), p(encountered), frag, omit, ing]
    zooEncFracs = [reshZooEncStats(:,1:2),zooFrag,zooOmit,zooIng];
    
    % Initialise output array on first iteration
    if iLoc == 1
        outZooEncStats = zeros([size(zooEncFracs),config.nLocs]);
    end
    outZooEncStats(:,:,iLoc) = zooEncFracs;
end
%%
for iFigure = 1:2
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.65 0.25],'Color','w')
    haxis = zeros(1,config.nLocs);

    for iLoc = 1:config.nLocs
        haxis(iLoc) = subaxis(1,config.nLocs,iLoc,'Spacing',0.035,'Padding',0.001,'Margin',0.11);
        ax(iLoc).pos = get(haxis(iLoc),'Position');
        ax(iLoc).pos(1) = ax(iLoc).pos(1) - 0.075;
        set(haxis(iLoc),'Position',ax(iLoc).pos)

        if iFigure == 1
            y = squeeze(outZooEncStats(:,3:5,iLoc));
            nStats = 3;
        elseif iFigure == 2
            y = squeeze(outZooEncStats(:,1:2,iLoc));
            nStats = 2;    
        end
    
        % Colour plot
        coloursStats = brewermap(nStats,'Paired');
        harea = area(haxis(iLoc),auxTermsDepths(:,iLoc),y);
        for iLt = 1:nStats
            harea(iLt).FaceColor = coloursStats(iLt,:);
        end
        
        ylim([0 1])
        set(haxis(iLoc), 'XDir', 'reverse');
        set(haxis(iLoc), 'XTick', [0 500 1000 1500 2000 3000 4000]);
        xtickformat('%d')
        if iLoc == 1
            xlabel('Depth (m)');
        end
        box on
        title(desiredModLocationOrder(iLoc),'FontSize',14);
        view([90 -90])
        
    end % iLoc

    % Legend
    if iFigure == 1
        lg = legend(harea,{'Fragmented','Omitted','Ingested'});
        lg.Position(2) = 0.75;
        % lg.Title.String = 'Mesozooplankton encounter outcome';
    elseif iFigure == 2
        lg = legend(harea,{'Evaluated','Encountered'});
        lg.Position(2) = 0.79;
        % lg.Title.String = 'Mesozooplankton grazing outcome';
    end
    lg.Position(1) = 0.82;
    lg.ItemTokenSize = [11,5];
    lg.FontSize = 10.5;
    set(lg,'Box','off');

    if iFigure == 1
        saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_zoostats_ingest'))
    elseif iFigure == 2
        saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_zoostats_encounters'))
    end

end % iFigure

end % plotEncounterTerms

% *************************************************************************

function plotRelativeContributionParticleTypes(figureSubfolderName,testrunDir,...
    figurePrefixName,config,filenameTimeseriesInformation,filenameObsPocFlux,...
    fullpathModelRunsDir,filenameSlamsOutput)

%% SECTION 0: Metadata

load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};

% Particle types
typeLabels = {'Faecal pellet',...
              'Aggregate',...
              'Living phyto. cell',...
              'Dead phyto. cell',...
              'TEP',...
              'Clay',...
              'Zoo. carcass',...
              'Inorganic fragment',...
              'Organic fragment'};
% typeLabels = {'Faecal pellet',...
%               'Aggregate',...
%               'Living phyto. cell',...
%               'Dead phyto. cell',...
%               'TEP',...
%               'Clay',...
%               'Zoo. carcass',...
%               'Inorganic fragment'};
nPartTypes = length(typeLabels);

% Colours
coloursTy = brewermap(nPartTypes,'Paired');

%% SECTION 1: Load modelled data

load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')

% Extrac particle type data (nTypes x nDephts x nAtts x nLocs)
partTypeData = output.particle.avgAttInTyAnnual;

%% SECTION 2: Transform into fractions

% Transform quantities into fractions
fracPartType = zeros(size(partTypeData,1,2,4));
for iLoc = 1:config.nLocs  
    for iZ = 1:size(config.availImagSysDeployDepths,1)
        partTypeNum = squeeze(partTypeData(:,iZ,1,iLoc)); % 1 = pnum
        partTypeDiam = squeeze(partTypeData(:,iZ,2,iLoc)); % 2 = diam
        totPartNum = sum(partTypeNum,'omitnan');
        for iTy = 1:size(partTypeData,1)
            fracPartType(iTy,iZ,iLoc) = partTypeNum(iTy)/totPartNum;
        end
    end
end

% Rearrange model locations
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'};
[~,reorderModLocIdx] = ismember(desiredLocationOrder,currentModLocationOrder); % get reordering indices
fracPartType = fracPartType(:,:,reorderModLocIdx);

%% SECTION 3: Depth selection for plotting

% POC flux compilation depths
load(fullfile('.','data','processed',filenameObsPocFlux),'LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED')

% Rearrange observed locations depths
[~,reorderObsLocIdx] = ismember(desiredLocationOrder,STATION_NAMES); % get reordering indices
LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED(:,reorderObsLocIdx,:,:);

% Calculate annual average
obsDepthsAnnualAvg = squeeze(mean(LOC_DEPTH_HORIZONS_POC_OBS_ADAPT_FILLED,1,'omitnan'));

% Replace the last obserevd depth in sediment traps with a depth more
% relevant for particle number observations
localImagedDepths = output.deployment.imagSys;
for iLoc = 1:config.nLocs
    if all(isnan(localImagedDepths(:,iLoc)))
        localImagedDepths(:,iLoc) = localImagedDepths(:,iLoc-1);
    end
    idxLastObsDepth = find(~isnan(localImagedDepths(:,iLoc)), 1, 'last'); % last non-NaN
    obsDepthsAnnualAvg(iLoc,:,4) = localImagedDepths(idxLastObsDepth,iLoc);
end

%% SECTION 4: Plot pie plot

% Extract modeled data at depths of interest
nImagedDepths = size(fracPartType,2);
nTargetDepths = size(obsDepthsAnnualAvg,3);
fracPartTypeDh = NaN(nPartTypes,nTargetDepths,config.nLocs); 
fracPartTypeDhDepths = NaN(nPartTypes,nTargetDepths,config.nLocs);

validDepthMask = false(nPartTypes,nImagedDepths,config.nLocs,nTargetDepths);
for iDh = 1:nTargetDepths
    for iLoc = 1:config.nLocs
        depthRangeTop = obsDepthsAnnualAvg(iLoc,1,iDh);    % upper boundary depth
        depthRangeBottom = obsDepthsAnnualAvg(iLoc,2,iDh); % lower boundary depth

        % Only proceed if both are numbers rather than NaN
        if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
            % Find the closest depth greater than or equal to depthRangeTop
            [~,idxTop] = min(abs(localImagedDepths(:,iLoc) - depthRangeTop));

            % Find the closest depth less than or equal to depthRangeBottom
            [~,idxBottom] = min(abs(localImagedDepths(:,iLoc) - depthRangeBottom));
            if idxBottom > nImagedDepths
                idxBottom = nImagedDepths;
            end

            % Extract valid modelled depth range
            validDepthMask(:,(idxTop:idxBottom),iLoc,iDh) = true;

            % Calculate the average monthly unique modelled depth
            fracPartTypeDhDepths(:,iDh,iLoc) =...
                (localImagedDepths(idxTop,iLoc)+localImagedDepths(idxBottom,iLoc))/2;
        end
    end % iLoc
end % iDh

% Extract particle type data
for iDh = 1:nTargetDepths
    currentMask = squeeze(validDepthMask(:,:,:,iDh));
    partTypeValues = fracPartType(:,:,:); 
    partTypeValues(~currentMask) = NaN;
    fracPartTypeDh(:,iDh,:) = squeeze(mean(partTypeValues,2,'omitnan'));
end

% Tranform into percentages of 100
fracPartTypeDhPct = 100.*fracPartTypeDh;

% Values for inset and main pie
threshold = 1; % 1%
mainTarget = 100 - threshold; % main pie target

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.60 0.60],'Color','w')
haxis = zeros(nTargetDepths,config.nLocs);

for iLoc = 1:config.nLocs
    for iDh = 1:nTargetDepths
        valsPct = squeeze(fracPartTypeDhPct(:,iDh,iLoc));

        % Only continue if there are data
        if ~isnan(sum(valsPct))
            [valsSorted, sortIdx] = sort(valsPct, 'descend');
            labelsSorted = typeLabels(sortIdx);
            
            % Accumulate until >= 99%
            cumVals = cumsum(valsSorted);
            cutIdx = find(cumVals >= mainTarget, 1, 'first');
            
            % Main pie values
            valsMain = valsSorted(1:cutIdx);
            labelsMain = labelsSorted(1:cutIdx);

            % If the last value overflows the target, split it
            overflow = cumVals(cutIdx) - mainTarget;
            if overflow > 0
                % Reduce last slice in main
                valsMain(end) = valsMain(end) - overflow;

                % Leftover goes into inset
                leftover = overflow;
                valsInset = [leftover; valsSorted(cutIdx+1:end)];
                labelsInset = [labelsMain(end)'; labelsSorted(cutIdx+1:end)'];

                % Reconfigure main labels
                valsMain = [valsMain; threshold];
                labelsMain = [labelsSorted(1:cutIdx)'; 'Other'];
            else
                valsInset = valsSorted(cutIdx+1:end);
                labelsInset = labelsSorted(cutIdx+1:end)';
            end
            
            % Normalise inset to 100%
            if ~isempty(valsInset)
                valsInsetRel = 100 * (valsInset / sum(valsInset));
            else
                valsInsetRel = [];
            end

            % Subplot position
            subplotIdx = (iDh-1)*config.nLocs + iLoc;
            haxis(iDh,iLoc) = subaxis(nTargetDepths,config.nLocs,subplotIdx,...
                'Spacing',0.005,'Padding',0.001,'Margin',0.05);
            
            ax(iDh,iLoc).pos = get(haxis(iDh,iLoc),'Position');
            ax(iDh,iLoc).pos(1) = ax(iDh,iLoc).pos(1) - 0.06;
            set(haxis(iDh,iLoc),'Position',ax(iDh,iLoc).pos)

            %%%%%%%%%%%%%%% Main pie %%%%%%%%%%%%%%%
            hMain = pie(haxis(iDh,iLoc), valsMain, repmat({''}, size(valsMain)));
            
            % Colour patches (every other handle is text)
            patchHandles = findobj(hMain, 'Type', 'Patch');
            
            % Assign colors based on typeLabels
            for k = 1:numel(patchHandles)
                % Find the index of this label in typeLabels
                lbl = labelsMain{k};
                if strcmp(lbl, 'Other')
                    % Gray color for "Other"
                    patchHandles(k).FaceColor = [1 1 1]; % white
                else
                    idx = find(strcmp(typeLabels, lbl));
                    patchHandles(k).FaceColor = coloursTy(idx,:);
                end
            end
            axis(haxis(iDh,iLoc),'equal');

            % Title
            if iDh == 1
                t = title(desiredLocationOrder{iLoc},'FontSize',16);
                pos = get(t,'Position'); pos(2) = pos(2)+0.12; set(t,'Position',pos);
            end

            %%%%%%%%%%%%%%% Inset pie for leftover small fractions %%%%%%%%%%%%%%%
            if sum(valsInsetRel) > 0
                insetPos = get(haxis(iDh,iLoc),'Position');
                insetPos = [insetPos(1)+0.05*insetPos(3), ... % move left
                            insetPos(2)+0.65*insetPos(4), ...
                            0.4*insetPos(3), 0.4*insetPos(4)];
                ax2 = axes('Position',insetPos);

                % Plot inset pie
                hInset = pie(ax2, valsInsetRel, repmat({''}, size(valsInsetRel)));

                patchInset = findobj(hInset,'Type','Patch');
                for k = 1:numel(patchInset)
                    lbl = labelsInset{k};
                    if strcmp(lbl, 'Other')
                        patchInset(k).FaceColor = [1 1 1]; % white
                    else
                        idx = find(strcmp(typeLabels,lbl));          % index in master type list
                        patchInset(k).FaceColor = coloursTy(idx,:);  % color from brewer map
                    end
                end
                axis(ax2,'equal');
            end
        end 
    end % iDh
end % iLoc

% Legend
maxLegendTypes = 9;  
nLegendTypes = min(maxLegendTypes, numel(typeLabels)); % cap at 9 or less

hLegend = gobjects(nLegendTypes,1); 
for i = 1:nLegendTypes
    hLegend(i) = patch(NaN,NaN,coloursTy(i,:),'DisplayName',typeLabels{i});
end

lgd = legend(hLegend, 'Location', 'eastoutside');
lgd.Title.String = 'Particle types';
lgd.FontSize = 11;
lgd.Position(1) = 0.88; 
lgd.Position(2) = 0.76;
lgd.ItemTokenSize = [11, 5]; 
lgd.Box = 'off';

saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_pie'))

%% SECTION 5: Plot histogram plot

% Desired depths to reveal
desiredDepths = [50,100,200,500,1000,1500,2000];
maxDesiredDepth = desiredDepths(end);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.60 0.25],'Color','w')
haxis = zeros(1,config.nLocs);

for iLoc = 1:config.nLocs
    haxis(1,iLoc) = subaxis(1,config.nLocs,iLoc,'Spacing',0.03,'Padding',0.001,'Margin',0.08);
    ax(iLoc).pos = get(haxis(iLoc),'Position');
    ax(iLoc).pos(1) = ax(iLoc).pos(1) - 0.038;
    set(haxis(iLoc),'Position',ax(iLoc).pos)

    % Extract data
    localFracPartType = squeeze(fracPartType(:,:,iLoc));

    % Discard the last sampled depth (seafloor)
    nanCols = all(isnan(localFracPartType), 1); % logical 1×58, true where all values in the column are NaN
    firstNanCol = find(nanCols, 1); % find the first such column
    if ~isempty(firstNanCol) && firstNanCol > 1
        localFracPartType(:,firstNanCol-1) = NaN; 
    end

    % Limit to maximum desired depth
    [~,idxMaxDepth] = min(abs(config.availImagSysDeployDepths - maxDesiredDepth));
    actualPnumDepths = config.availImagSysDeployDepths(1:idxMaxDepth);

    % Generate normalised depth axis
    normPnumDepth = linspace(1, idxMaxDepth, idxMaxDepth);

    % Reshape for plotting
    reshLocalFracType = 100.*localFracPartType(:,1:idxMaxDepth)'; % transpose localFracPartType from nTypes x nDepths --> nDepths x nTypes
    reshLocalFracType(reshLocalFracType==0) = NaN;
    
    % Plot
    hbarh = barh(haxis(iLoc),normPnumDepth,reshLocalFracType,'stacked'); 
    hold on
    for iTy = 1:nPartTypes
        hbarh(iTy).FaceColor = coloursTy(iTy,:);
    end
    xlim([0 100])
    set(haxis(iLoc), 'YDir', 'reverse');

    hold off

    % Adjust y-axis limits to remove extra space above and below the bars
    ylim(haxis(iLoc), [min(normPnumDepth)-0.5, max(normPnumDepth)+0.5]);

    % Calculate y-tick positions based on closest depths
    [~,idxDesiredDepths] = arrayfun(@(d) min(abs(actualPnumDepths - d)), desiredDepths);
    ytickPositions = normPnumDepth(idxDesiredDepths);
    yticks(haxis(iLoc), ytickPositions);
    yticklabels(haxis(iLoc), arrayfun(@num2str, desiredDepths, 'UniformOutput', false));

    if iLoc == 1
        ylabel('Depth (m)');
    end

    % xlabel('Particle types')
    title(desiredLocationOrder(iLoc),'FontSize',12);

end % iLoc

% Legend
lg = legend(hbarh,typeLabels);
lg.Position(1) = 0.89; 
lg.Position(2) = 0.54;
lg.ItemTokenSize = [11,5];
lg.FontSize = 11;
set(lg,'Box','on');

saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_hist'))
   
end % plotRelativeContributionParticleTypes

% *************************************************************************

function plotLocalForcingDataDepthDistributed(dailyDepthDistribData,...
    dailyDepthDistribDataMetadata,config,monthStart,desiredModLocationOrder,...
    stationTags,figureSubfolder,testrunDir)

% Mapping properties
myColourmap = brewermap(1000,'*Spectral'); % Spectral is the equivalent of Jet
x = 1:365;

for iLoc = 1:config.nLocs

    locDepths = config.gridDepths(:,iLoc);
    validDepthIdx = find(~isnan(locDepths));
    validDepths   = locDepths(validDepthIdx);

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.35 0.55],'Color','w')
    haxis = zeros(size(dailyDepthDistribDataMetadata,1),1);

    for iDataset = 1:size(dailyDepthDistribDataMetadata,1)
        haxis(iDataset) = subaxis(3,3,iDataset,'Spacing',0.025,'Padding',0.005,'Margin',0.10);
        ax(iDataset).pos = get(haxis(iDataset),'Position');
        if (iDataset > 3 && iDataset < 7)
            ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.025;
        end
        if (iDataset >=7)
            ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.050;
        end
        set(haxis(iDataset),'Position',ax(iDataset).pos) 

        datasetName = dailyDepthDistribDataMetadata{iDataset,1};
        titleStr = dailyDepthDistribDataMetadata{iDataset,2}; 
        caxisMax = dailyDepthDistribDataMetadata{iDataset,4}; 
        caxisMin = dailyDepthDistribDataMetadata{iDataset,3}; 
        data = dailyDepthDistribData.(datasetName); 

        y = validDepths;
        dataMatrix = squeeze(data(validDepthIdx,:,iLoc));

        % Plot
        h = imagesc(x,y,dataMatrix);
        % set(gca, 'YScale', 'log'); % Set y-axis to log scale
        set(h, 'AlphaData', ~isnan(dataMatrix));  % makes NaNs fully transparent
        colormap(myColourmap)
        caxis([caxisMin, caxisMax]);
        shading flat; % smooth the color transitions
        hold on

        % xticks
        xticks(monthStart);
        xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
        xlim([1 365]);
        axh = gca;
        axh.XAxis.FontSize = 8;
        set(gca, 'XTickLabelRotation', 0); % keep x-tick labels horizontal

        % yticks
        ytickValues = [100, 500, 1000, 1500, 2000];
        ytickValues = ytickValues(ytickValues <= max(validDepths));
        yticks(ytickValues); 

        % Show labels only for selected months
        if ismember(iDataset, [1, 4, 7])
            yticklabels(arrayfun(@num2str, ytickValues, 'UniformOutput', false));
        else
            yticklabels([]); % hide y-axis labels
        end
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.02, 0.02]; % make tick marks longer

        % Title
        title(titleStr,'FontSize',12,'FontWeight', 'normal')

        cb = colorbar('Location','eastoutside');
        % cb.Ticks = ticks;
        % cb.TickLabels = tickLabels;
        cb.FontSize = 11; % font size of tick labels
        % cb.Label.FontSize = 18; % font size of colour bar string
        % cb.Position(1) = cb.Position(1) + 0.09;
        % cb.Position(2) = cb.Position(2) - 0.60;
        % cb.Position(3) = 0.025; % WIDTH
        % cb.Position(4) = 0.74; % LENGTH
        % cb.Label.String = titleStr;

    end % iDataset

    a = axes; % create a new axis
    t = title(desiredModLocationOrder(iLoc),'FontSize',16);
    yl = ylabel('Depth (m)','FontSize',11);
    a.Visible = 'off';
    t.Visible = 'on';
    t.Position(1) = t.Position(1) - 0.07; t.Position(2) = t.Position(2) + 0.035;
    yl.Visible = 'on';
    yl.Position(1) = yl.Position(1) - 0.06;
    yl.Position(2) = yl.Position(2) - 0.05;
      
    % Save figure 
    saveFigureInFolder(figureSubfolder,strcat(testrunDir,'_localforcing_depthdistrib_',stationTags{iLoc}))

end % iLoc

end % plotLocalForcingDataDepthDistributed

% *************************************************************************

function plotLocalForcingDataSurfaceDistribution(dailySurfaceData,dailySurfaceDataMetadata,...
    config,monthStart,desiredModLocationOrder,stationTags,figureSubfolder,testrunDir)

colourFirstFour = brewermap(4,'*Paired');
shiftRight = 0.08;
heigtSubplot = 0.28;
x = 1:365;

for iLoc = 1:config.nLocs
        
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.30],'Color','w') 
    
    for iSubplot = 1:3
        
        % NPP and chla
        if (iSubplot == 1) 
    
            % Extract NPP data
            datasetName = dailySurfaceDataMetadata{1,1};
            nppLabelStr = dailySurfaceDataMetadata{1,2};
            nppdata = dailySurfaceData.(datasetName); 

            % Extract chla data
            datasetName = dailySurfaceDataMetadata{2,1};
            chlaLabelStr = dailySurfaceDataMetadata{2,2};
            chladata = dailySurfaceData.(datasetName); 
    
            ax1 = axes('Position', [shiftRight 0.65 0.80 heigtSubplot]);
    
            yyaxis left
            lineColorL = colourFirstFour(1,:); 
            plot(ax1,x,nppdata(:,iLoc),'-','Color',lineColorL,'LineWidth',2); hold on;
            if iLoc < 6
                set(ax1,'YScale','log');
                set(ax1,'YTick', [200, 500, 1000, 1500]);
                ylim([200 2000])
            else
                set(ax1,'YTick', [50, 100, 200]);
                ylim([20 250])
            end
            xticks(monthStart);
            xticklabels([])
            xlim([1 365]);
            ylabel(nppLabelStr,'FontSize',11,'Color','k')
            ax1.TickDir = 'in';
            ax1.YColor = lineColorL; 
            ax1.XGrid = 'on';
            ax1.YGrid = 'off'; 
    
            yyaxis right
            lineColorR = colourFirstFour(2,:); 
            plot(ax1,x,chladata(:,iLoc),'-','Color',lineColorR,'LineWidth',2); hold on;
            if iLoc < 6
                set(ax1,'YScale','log');
                set(ax1,'YTick', [0.01, 0.1, 0.2, 0.5, 1]);
                ylim([0.01 1])
            else
                set(ax1,'YScale','log');
                set(ax1,'YTick', [0.1, 0.2, 1, 2]);
                ylim([0.1 2.1])
            end
            ytickformat('%.1f');
            xticks(monthStart);
            xticklabels([])
            xlim([1 365]);
            ylabel(chlaLabelStr,'FontSize',11,'Color','k')
            ax1.TickDir = 'in';
            ax1.YColor = lineColorR; 
            ax1.XGrid = 'on';
            ax1.YGrid = 'off'; 
    
            box on
    
        % MLD and PAR0
        elseif (iSubplot == 2) % C products
    
            % Extract PAR0 data
            datasetName = dailySurfaceDataMetadata{3,1};
            par0LabelStr = dailySurfaceDataMetadata{3,2};
            par0data = dailySurfaceData.(datasetName); 

            % Extract MLD data
            datasetName = dailySurfaceDataMetadata{4,1};
            mldLabelStr = dailySurfaceDataMetadata{4,2};
            mlddata = dailySurfaceData.(datasetName); 
              
            ax2 = axes('Position', [shiftRight 0.30 0.80 heigtSubplot]);
    
            yyaxis left
            lineColorL = colourFirstFour(3,:); 
            plot(ax2,x,par0data(:,iLoc),'-','Color',lineColorL,'LineWidth',2); hold on;
            ylim([0 180])
            ylabel(par0LabelStr,'FontSize',11,'Color','k')
            xticks(monthStart);
            xticklabels([])
            xlim([1 365]);
            ax2.TickDir = 'in';
            ax2.XGrid = 'on';
            ax2.YGrid = 'off'; 
            ax2.YColor = lineColorL; % Set the color of the left y-axis to black
            
            yyaxis right
            lineColorR = colourFirstFour(4,:); 
            plot(ax2,x,mlddata(:,iLoc),'-','Color',lineColorR,'LineWidth',2); hold on;
            if iLoc < 6
                ylim([0 320])
            elseif iLoc == 6
                ylim([0 550])
            end
            ylabel(mldLabelStr,'FontSize',11,'Color','k')
            xticks(monthStart);
            xticklabels([])
            xlim([1 365]);
            ax2.TickDir = 'in';
            ax2.XAxis.FontSize = 10.5;
            ax2.XGrid = 'on';
            ax2.YGrid = 'off'; 
            ax2.YColor = lineColorR;
    
            box on
            
        % Dust
        elseif (iSubplot == 3) 
    
            % Extract data
            datasetName = dailySurfaceDataMetadata{5,1};
            dustLabelStr = dailySurfaceDataMetadata{5,2};
            dustdata = dailySurfaceData.(datasetName);
    
            ax3 = axes('Position', [shiftRight 0.10 0.80 0.14]); 
           
            plot(ax3,x,dustdata(:,iLoc),'-','Color','r','LineWidth',2); hold on;
            set(ax3,'YScale','log');
            set(ax3,'YTick', [1e-3, 1e-2, 1e-1, 1]);
            ylim([1e-3 7])
            ylabel(dustLabelStr)
            xticks(monthStart);
            xticklabels({'J','F','M','A','M','J','J','A','S','O','N','D'})
            xlim([1 365]);
            ytickformat('%.1f');
            ax3.XGrid = 'on';
            ax3.YGrid = 'off'; 
            hold on
    
            box on
    
        end
       
    end

    title(ax1, desiredModLocationOrder{iLoc}, 'FontSize', 16);

    % Save figure 
    saveFigureInFolder(figureSubfolder,strcat(testrunDir,'_localforcing_surfacedistrib_',stationTags{iLoc}))
      
end % iLoc

end % plotLocalForcingDataSurfaceDistribution

% *************************************************************************

function plotPftProbabilityDistributions(figureSubfolderName,testrunDir,...
    figurePrefixName,config,fullpathModelRunsDir,filenameSlamsOutput)

%% SECTION 1: Load data

load(fullfile(fullpathModelRunsDir,filenameSlamsOutput),'output')

% Extract monthly aux. data (4 x 12 x nLocs)
monthlyPftBiomass = squeeze(output.aux.monthly(config.aux.idxDiatBiomass:config.aux.idxPicoBiomass,1,:,:));  
monthlyPftRelatProb = squeeze(output.aux.monthly(config.aux.idxAvgProbDiat:config.aux.idxAvgProbPico,1,:,:)); 

% Extract seasonal aux. data (4 x nDepths x 4 x nLocs)
seasonalPftProbDepthDistrib = squeeze(output.aux.seasonal(config.aux.idxProfProbDiat:config.aux.idxProfProbPico,:,:,:));
 
% Rearrange model locations
currentModLocationOrder = {'EqPac','HOT/ALOHA','BATS/OFP','PAP-SO','OSP','HAUSGARTEN'};
desiredModLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
[~,reorderModLocIdx] = ismember(desiredModLocationOrder,currentModLocationOrder); % get reordering indices

% Reorder to match desired locations
monthlyPftBiomass   = monthlyPftBiomass(:,:,reorderModLocIdx);
monthlyPftRelatProb = monthlyPftRelatProb(:,:,reorderModLocIdx);
seasonalPftProbDepthDistrib = seasonalPftProbDepthDistrib(:,:,:,reorderModLocIdx);
seasonalPftProbDepthDistrib(seasonalPftProbDepthDistrib==0) = NaN;

% Get depths to plot seasonal data
validMask = seasonalPftProbDepthDistrib > 0 & isfinite(seasonalPftProbDepthDistrib);
validDepth = squeeze(any(any(any(validMask, 1), 3), 4)); % collapse across PFTs, seasons, and locations
maxDepthIdx = find(validDepth, 1, 'last'); % deepest depth index with any data
auxTermsDepths = config.gridDepths(1:maxDepthIdx,reorderModLocIdx);

% Recalc to get relative biomass abundance
totalBiomass = sum(monthlyPftBiomass, 1);
totalBiomass(totalBiomass == 0) = NaN;
monthlyRelatPftBiomass = monthlyPftBiomass ./ totalBiomass;
monthlyRelatPftBiomass(isnan(monthlyRelatPftBiomass)) = 0;

seasonLabel = {'Winter','Spring','Summer','Autumn'};

% PFTs
typeLabels = {'Diatom',...
              'Large flagellate',...
              'Coccolithophore',...
              'Picophytoplankton'};
nPfts = length(typeLabels);

% Colours
% coloursTy = brewermap(nPfts,'Accent');
% coloursTy = coloursTy([3 2 4 1], :); % reorder
hexColors = {'#F9B232', '#268874', '#9D9D9C', '#9EC043'};
hex2rgb = @(hex) sscanf(hex(2:end),'%2x%2x%2x',[1 3]) / 255;
rgbColors = cellfun(hex2rgb, hexColors, 'UniformOutput', false);
coloursTy = vertcat(rgbColors{:});

%% SECTION 2: Plot pie plot summaries

monthNames = config.mappingProps.labelMonths;

for iFigure = 1:2
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.75],'Color','w')
    haxis = zeros(config.nLocs,12);

    for iLoc = 1:config.nLocs
        for iMonth = 1:12
            if (iFigure == 1)
                vals = squeeze(monthlyPftRelatProb(:,iMonth,iLoc));
                gentitle = 'Normalised seeding probability (avg euphotic layer)';
                figureLabel = strcat(figurePrefixName,'_seedingprob_pie');
            elseif (iFigure == 2)
                vals = squeeze(monthlyRelatPftBiomass(:,iMonth,iLoc));
                gentitle = 'Relative biomass abundance';
                figureLabel = strcat(figurePrefixName,'_relatbiom_pie');
            end 
    
            % Only continue if there are data
            if ~isnan(sum(vals))
                [valsSorted, sortIdx] = sort(vals, 'descend');
                labelsSorted = typeLabels(sortIdx);
    
                % Main pie values
                valsMain = valsSorted;
                labelsMain = labelsSorted;
    
                % Subplot position
                subplotIdx = (iMonth-1)*config.nLocs + iLoc;
                haxis(iMonth,iLoc) = subaxis(12,config.nLocs,subplotIdx,...
                    'Spacing',0.005,'Padding',0.001,'Margin',0.07);
                
                ax(iMonth,iLoc).pos = get(haxis(iMonth,iLoc),'Position');
                ax(iMonth,iLoc).pos(1) = ax(iMonth,iLoc).pos(1) - 0.03;
                set(haxis(iMonth,iLoc),'Position',ax(iMonth,iLoc).pos)
    
                if sum(valsMain) > 0
                    hMain = pie(haxis(iMonth,iLoc), valsMain, repmat({''}, size(valsMain)));
                
                    % Colour patches
                    patchHandles = findobj(hMain, 'Type', 'Patch');
                    for k = 1:numel(patchHandles)
                        lbl = labelsMain{k};
                        idx = find(strcmp(typeLabels, lbl));
                        patchHandles(k).FaceColor = coloursTy(idx,:);
                    end
                else
                    % Nothing to plot – keep empty space
                    axis(haxis(iMonth,iLoc),'off');
                end

                % Title
                if iMonth == 1
                    t = title(desiredModLocationOrder{iLoc},'FontSize',13);
                    pos = get(t,'Position'); pos(2) = pos(2)+0.12; set(t,'Position',pos);
                end
                % Month label on the left of the first column
                if iLoc == 1
                   text(-0.15, 0.5, monthNames{iMonth}, ...
                        'Units','normalized', ...
                        'FontSize',11, ...
                        'HorizontalAlignment','right', ...
                        'VerticalAlignment','middle'); 
                end
            end 

        end % iMonth
    end % iLoc

    % Legend
    hLegend = gobjects(nPfts,1);
    for iPft = 1:nPfts
        hLegend(iPft) = patch(NaN,NaN,coloursTy(iPft,:),...
            'DisplayName',typeLabels{iPft});
    end
    lgd = legend(hLegend,typeLabels,'Orientation','horizontal','FontSize',11);
    lgd.Box = 'off';
    lgd.ItemTokenSize = [15, 25]; 
    lgd.Position = [0.15 0.01 0.70 0.05];

    % --- Figure-level title (version-independent) ---
    a = axes; % create a new axis
    a.Visible = 'off';
    t = title(gentitle,'FontSize',16,'FontWeight','bold');
    t.Visible = 'on';
    t.Position(1) = t.Position(1) - 0.07; t.Position(2) = t.Position(2) + 0.050;

    saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figureLabel))
end

%% SECTION 3: Plot depth distribution of relative seeding probabilities of PFTs

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.50],'Color','w')
haxis = zeros(config.nLocs*4,1);

for iSeason = 1:4
    for iLoc = 1:config.nLocs

        % Compute linear subplot index (row-wise fill)
        iSubplot = (iSeason - 1) * config.nLocs + iLoc;
        haxis(iSubplot) = subaxis(4,config.nLocs,iSubplot,'Spacing',0.038,'Padding',0,'Margin',0.14);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
  
        %(4 PFTs x nDepths x 4 x nLocs)
        for iPft = 1:nPfts
            plot(haxis(iSubplot), ...
                squeeze(seasonalPftProbDepthDistrib(iPft,1:maxDepthIdx,iSeason,iLoc)), ...
                auxTermsDepths, ...
                'LineStyle', '-', ...
                'Color', coloursTy(iPft,:), ...
                'LineWidth', 2); hold on;
        end
        hold off;

        xlim([0 0.7])
        xTickValues = 0:0.30:0.70;
        xticks(xTickValues)
        xticklabels(xTickValues)

        ylim([0 150])
        yticks([0 50 100 150])
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.03, 0.03]; % make tick marks longer

        if (iLoc == 1)
            yticklabels({'0','50','100','150'})
        else
            yticklabels([])
        end

        set(gca,'YDir','Reverse','XAxisLocation','Bottom','xlabel',[],'ylabel',[],'FontSize', 10)

        if iSeason == 1
            text(0.5, 1.35, desiredModLocationOrder{iLoc}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 14, 'FontWeight', 'bold');
        end

        % Shift all plots a little bit up
        if (iSubplot <= 6)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.040; 
        elseif (iSubplot > 6 && iSubplot <= 12)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)+0.015; 
        elseif (iSubplot > 12 && iSubplot <= 18)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.010; 
        elseif (iSubplot > 18)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2)-0.035; 
        end
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.065; 
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

    end % iLoc
end % iSeason

hLegend = gobjects(nPfts,1);
for iPft = 1:nPfts
    hLegend(iPft) = patch(NaN,NaN,coloursTy(iPft,:),...
        'DisplayName',typeLabels{iPft});
end
lgd = legend(hLegend,typeLabels,'Orientation','vertical','FontSize',11);
lgd.Position(1) = 0.80; 
lgd.Position(2) = 0.11;
lgd.ItemTokenSize = [11,5];
lgd.FontSize = 11;
set(lgd,'Box','on');

% Add season labels above each row
for iSeason = 1:4
    % Get subplot indices for the current row
    rowIdx = (iSeason - 1) * config.nLocs + (1:config.nLocs);

    % Extract the axis positions for the row
    rowPositions = arrayfun(@(h) get(haxis(h), 'Position'), rowIdx, 'UniformOutput', false);
    rowPositions = cat(1, rowPositions{:});  % Convert to numeric matrix: 6x4 [x y w h]

    % Compute average x center of the row (between subplots 3 and 4 is safe)
    xCenters = rowPositions(:,1) + rowPositions(:,3)/2;
    xCenterRow = mean(xCenters);

    % Compute highest y position in that row (top of the tallest axis)
    topEdges = rowPositions(:,2) + rowPositions(:,4);
    yTopRow = max(topEdges) + 0.008;  % Add small offset above the top

    % Add season label annotation
    annotation('textbox', [xCenterRow - 0.05, yTopRow, 0.1, 0.03], ...
        'String', seasonLabel{iSeason}, ...
        'FontWeight', 'bold', ...
        'FontSize', 12, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center');
end

% Give common xlabel, ylabel and title to your figure
% Create a new axis
a = axes;
xl = xlabel('Relative PFT seeding probability','FontSize',14);
yl = ylabel('Depth (m)','FontSize',14);

% Specify visibility of the current axis as 'off'
a.Visible = 'off';
% Specify visibility of Title, XLabel, and YLabel as 'on'
xl.Visible = 'on';
yl.Visible = 'on';
yl.Position(1) = yl.Position(1) - 0.09; 
xl.Position(1) = 0.38; xl.Position(2) = xl.Position(2) - 0.04;

saveFigureInFolder(figureSubfolderName,strcat(testrunDir,'_',figurePrefixName,'_depthdistrib'))

end % plotPftProbabilityDistributions

% *************************************************************************
