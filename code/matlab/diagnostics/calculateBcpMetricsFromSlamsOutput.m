% function calculateBcpMetricsFromSlamsOutput(fullpathModelRunsDir,...
%     fullpathModelInputDataDir, fullpathProcessedDataDir, filenameOutputSlams,...
%     filenameSlamsRunGrid, filenameSlamsNumDepthLayers, filenameInputNpp,...
%     filenameInputMask, filenameInputZeu, filenameOutputSlamsBcpMetrics,...
%     choiceTypeGridDomain)

% % CALCULATEBCPMETRICSFROMSLAMSOUTPUT Calculates BCP metrics (Teff, PEeff, 
% % b, z* and xi) for SLAMS output.                                           
% 
% % Folder paths:
  fullpathModelRunsDir      = './tests/LOCALTS6/modelruns/';
  fullpathModelInputDataDir = './tests/LOCALTS6/modelinputdata/';
  fullpathProcessedDataDir  = './data/processed/';
% 
% % Filenames:
  filenameOutputSlams           = 'runsoutput.mat';
  filenameSlamsRunGrid          = 'grid_run.mat';
  filenameSlamsNumDepthLayers   = 'waterColNumDepthLayers.txt';
  filenameInputNpp              = 'npp_bicep.mat';
  filenameInputMask             = 'mask_custom_icefrac_cmems_chla_occci.mat';
  filenameInputZeu              = 'zeu_calculated_chlaoccci_mldifremer_pointonepercentpar0.mat';
  filenameOutputSlamsBcpMetrics = 'bcpmetrics.mat';
  choiceTypeGridDomain          = 2; % 1=global, 2=local

addpath(genpath(fullfile('.','modelresources','internal'))) 
addpath(genpath(fullfile('.','modelresources','external'))) 
addpath(genpath(fullfile('.','code','matlab')))
addpath(genpath(fullfile('.','data','raw'))) 
addpath(genpath(fullfile('.','data','processed'))) 

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS: CONFIGURATION PARAMETERS
% -------------------------------------------------------------------------

% Choices to be applied to sediment trap/radionuclide data when fitting b and z* 
isFluxNormalised = 0; % 0=do not normalise POC flux values to value at zref, 1=do normalise
choiceZref = 2;       % reference depth, 1=closest value to 100, 2=zeu, 3=inflexion point

% Load configuration parameters used in the model
config = loadModelConfigurationParameters(fullpathModelInputDataDir,...
    filenameSlamsRunGrid,filenameSlamsNumDepthLayers,choiceTypeGridDomain);

% Log progress
logID = fopen(fullfile(fullpathModelRunsDir,'logCalculateBcpMetric.log'),'w');

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - EXTRACT AND ORGANISE MODEL OUTPUT DATA
% -------------------------------------------------------------------------

% Upload model output data
load(fullfile(fullpathModelRunsDir,filenameOutputSlams),'output')

% Extract
tmpOutput = extractModelOutput(output);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - GET KEY DEPTH LEVELS FOR PLOTTING GLOBAL OUTPUT
% -------------------------------------------------------------------------

% Extract supporting data for calculations
[seqZeuMonthly,seqZeuAnnual,seqNppMonthly] = extractOutputSupportingData(...
    config,filenameInputZeu,filenameInputNpp,filenameInputMask,...
    choiceTypeGridDomain,logID);

% Define key depths 
[keyDepthLevelsMonthlyData,keyDepthLevelsAnnualData] = ...
    defineKeyDepthsForPlottingOutput(config,seqZeuMonthly,choiceTypeGridDomain);

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
% SECTION 5 - EXTRACT DEPTH HORIZON DATA
% -------------------------------------------------------------------------

% Add fractional uncertainty to arrays and format them for further processing
seqFluxMonthlyProfile = prepareSeqData(tmpOutput,fluxSimulatedRelUncertaintyOverall,'flux');
seqPnumMonthlyProfile = prepareSeqData(tmpOutput,pNumSimulatedRelUncertaintyOverall,'particleNumbers');

% Extract flux at key depth horizons and propagate uncertainty (in flux units)
[seqFluxMonthlyDh,~,~,~] = extractKeyDepthHorizonDataInModelledFluxes(...
    config,seqFluxMonthlyProfile,keyDepthLevelsMonthlyData);

% Extract particle numbers at key depth horizons and propagate uncertainty (in particle number units)
[seqPnumMonthlyDh,~,~,~] = extractKeyDepthHorizonDataInModelledParticleNumbers(...
    config,seqPnumMonthlyProfile,keyDepthLevelsMonthlyData);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 6 - COMPUTE BCP METRICS
% -------------------------------------------------------------------------

isSlurm = ~isempty(getenv('SLURM_JOB_ID'));

if isSlurm
    % Start parallel pool once in the main script
    c = parcluster('local');
    nSlurm = str2double(getenv('SLURM_CPUS_PER_TASK'));
    if isnan(nSlurm)
        nSlurm = feature('numcores');
    end
    numWorkers = min(nSlurm, feature('numcores'));
    if isempty(gcp('nocreate'))
        parpool(c, numWorkers);
    end
    fprintf('Parallel pool started with %d workers.\n', numWorkers);
else
    fprintf('Running in serial mode.\n');
end

tic; % start timing the processing

% .........................................................................

% Calculate Teff and PEeff

% Extract and group data for average and relative error values
baseData = { % 1st col: avg, 2nd col: err
    squeeze(seqFluxMonthlyDh(1,:,:,1,1)), squeeze(seqFluxMonthlyDh(1,:,:,2,1)); % POC flux at zeu
    squeeze(seqFluxMonthlyDh(2,:,:,1,1)), squeeze(seqFluxMonthlyDh(2,:,:,2,1)); % POC flux at zmeso
    squeeze(seqNppMonthly(:,:,1)),        squeeze(seqNppMonthly(:,:,2))
};

groupData = cell(3,1);
for i = 1:3
    groupData{i} = NaN(12,config.nLocs,2);
    groupData{i}(:,:,1) = baseData{i,1}; % avg
    groupData{i}(:,:,2) = baseData{i,2}; % err
end

% Assign to specific variables for clarity
[arrayFluxZeu,arrayFluxZmeso,arrayNpp] = deal(groupData{:});

% Calculations
fprintf('\nInitiate calculation of Teff...\n');
[teffAnnual,teffMonthly] = calculateEfficiencyMetrics('Teff',arrayFluxZmeso,arrayFluxZeu);
fprintf('\n...done.\n');
fprintf('\nInitiate calculation of PEeff...\n');
[peeffAnnual,peeffMonthly] = calculateEfficiencyMetrics('PEeff',arrayFluxZeu,arrayNpp);
fprintf('\n...done.\n');

% Store results for Teff and PEeff using a struct update loop
effMetrics = {'teff','peeff'};
annualResults = {teffAnnual, peeffAnnual};
monthlyResults = {teffMonthly, peeffMonthly};
for i = 1:length(effMetrics)
    thisMetric = effMetrics{i};
    metrics.(thisMetric).monthly = monthlyResults{i}(:,:,:);
    metrics.(thisMetric).annual  = annualResults{i}(:,:);
end

% .........................................................................

% Fit b, z* and xi

% Prepare POC flux data
pocFlux = squeeze(seqFluxMonthlyProfile(:,:,:,1,1));
pocFluxErr = squeeze(seqFluxMonthlyProfile(:,:,:,2,1)).*squeeze(seqFluxMonthlyProfile(:,:,:,1,1)); % from fractional to value
sedTrapDeployDepths = permute(repmat(output.deployment.sedTrap, [1, 1, 12]), [1, 3, 2]); % from 2D to 3D
depthBoundaries = permute(keyDepthLevelsMonthlyData,[4,1,2,3]);

% Extract flux data from zeu to zmeso
arrayFlux   = NaN(numel(config.availSedTrapDeployDepths),12,config.nLocs,2); % 4th dim: 1=avg, 2=err 
arrayDepths = NaN(numel(config.availSedTrapDeployDepths),12,config.nLocs);
for iLoc = 1:config.nLocs
    for iMonth = 1:12
        [selectedDepths,selectedFluxes,selectedErrors] = ...
            extractDataFromZrefToZmeso(sedTrapDeployDepths(:,iMonth,iLoc),...
                pocFlux(:,iMonth,iLoc),pocFluxErr(:,iMonth,iLoc),...
                choiceZref,squeeze(depthBoundaries(:,:,iMonth,iLoc)),...
                config.maxZeu);
        nEntriesSelected = numel(selectedDepths);
        if (~isempty(nEntriesSelected) && nEntriesSelected > 0) 
            if isFluxNormalised
                selectedFluxes = selectedFluxes./selectedFluxes(1);
                selectedErrors = selectedErrors./selectedErrors(1);
            end
            arrayFlux(1:nEntriesSelected,iMonth,iLoc,1) = selectedFluxes; 
            arrayFlux(1:nEntriesSelected,iMonth,iLoc,2) = selectedErrors;
            arrayDepths(1:nEntriesSelected,iMonth,iLoc) = selectedDepths;
        end
    end
end

% Prepare particle number data
scalingFactor = 1e3 ./ reshape(config.particleDiameterClassWidth, [config.nDiameterClasses, 1, 1, 1]);
seqPnumMonthlyDh = seqPnumMonthlyDh .* scalingFactor; % # part. L-1 --> # part. m-3 um-1

% Extract particle number data from class with largest particle number to class with lowest
% particle number
arrayPnum = NaN(size(seqPnumMonthlyDh));
arraySizeClassMiddle = NaN(numel(config.particleDiameterClassMiddle),4,12,config.nLocs);
for iLoc = 1:config.nLocs
    for iMonth = 1:12
        for iDh = 1:4
            dataAvg     = seqPnumMonthlyDh(:,iDh,iMonth,iLoc,1); % # part. m-3 um-1
            dataErr     = seqPnumMonthlyDh(:,iDh,iMonth,iLoc,2); % # part. m-3 um-1
            dataClasses = config.particleDiameterClassMiddle; % um
    
            % Define y0, iLastClass, iFirstClass
            y0 = max(dataAvg(:)); % # part. m-3 um-1
            iFirstClass = find(dataAvg(:) == y0, 1, 'first');
            iLastClass = find(dataAvg(:) ~= 0, 1, 'last'); % the function goes from max to end
            nClasses = iLastClass - iFirstClass + 1;

            % Proceed if at least we have 3 bins with data
            if nClasses < 3
                continue; % skip to next iDh
            end

            % Trim data ranges
            arrayPnum(1:nClasses,iDh,iMonth,iLoc,1) = dataAvg(iFirstClass:iLastClass); % # part. m-3 um-1
            arrayPnum(1:nClasses,iDh,iMonth,iLoc,2) = dataErr(iFirstClass:iLastClass); % # part. m-3 um-1
            arraySizeClassMiddle(1:nClasses,iDh,iMonth,iLoc) = dataClasses(iFirstClass:iLastClass); % um
        end
    end
end

% Start calculations
fprintf('\nInitiate calculation of b...\n');
[martinbAnnual,martinbMonthly] = calculateSpectrumMetrics('b',arrayFlux,arrayDepths);
fprintf('\n...done.\n');
fprintf('\nInitiate calculation of z*...\n');
[zstarAnnual,zstarMonthly] = calculateSpectrumMetrics('zstar',arrayFlux,arrayDepths);
fprintf('\n...done.\n');
fprintf('\nInitiate calculation of xi...\n');
[xiAnnual,xiMonthly] = calculateSpectrumMetrics('xi',arrayPnum,arraySizeClassMiddle);
fprintf('\n...done.\n');

% Store results for b and z* using a struct update loop
specMetrics = {'martinb','zstar'};
annualResults = {martinbAnnual, zstarAnnual};
monthlyResults = {martinbMonthly, zstarMonthly};
for i = 1:length(specMetrics)
    thisMetric = specMetrics{i};
    metrics.(thisMetric).monthly = monthlyResults{i}(:,:,:);
    metrics.(thisMetric).annual  = annualResults{i}(:,:);
end

% Store results for xi
metrics.xi.monthly = xiMonthly(:,:,:);
metrics.xi.annual  = xiAnnual(:,:);

% .........................................................................

% Finalise and close files
vars = who('metrics');
save(fullfile(fullpathModelRunsDir,filenameOutputSlamsBcpMetrics),vars{:},'-v7.3') 
fprintf(logID, 'Processing completed in %.2f seconds.\n', toc);    
fclose(logID);

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

function [qZeuMonthly,qZeuAnnual,qNppMonthly] = extractOutputSupportingData(...
    config,filenameInputZeu,filenameInputNpp,filenameInputMask,...
    choiceTypeGridDomain,logID)

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
        if (choiceTypeGridDomain == 1) % global
            iLat = config.iyBb(iLoc);
            iLon = config.ixBb(iLoc);
        elseif (choiceTypeGridDomain == 2) % local
            iLat = iLoc;
            iLon = iLoc;
        end
        
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

% end % calculateBcpMetricsFromSlamsOutput
