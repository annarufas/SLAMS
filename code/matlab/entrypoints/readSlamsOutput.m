function readSlamsOutput(fullpathModelRunsDir, fullpathModelInputDataDir,...
    filenameOutputSlams, filenameInputRunGrid, filenameInputNumDepthLayers,...
    choiceTypeGridDomain, choiceTestModelEquilibrium, choiceSaveDiskSpace,...
    choiceSeasonalOutput)

% READSLAMSOUTPUT Reads SLAMS output .bin files, extracts the data, processes 
% the data to get seasonal and annual averages, and produces .mat arrays for 
% further processing (BCP metric calculation) and plotting. Monthly values 
% are extracted for the last model year run (at equilibrium). Seasonal and 
% annual averages are calculated from weighted monthly values. The output
% arrays generated are:
%
% - Fluxes (mg m-2 d-2): POC, TEP-C, CaCO3, bSi and clay at sediment trap 
%   deployment depths and the seafloor.
% - Particle attributes (sorted by size class, velocity class, particle
%   type)
%      - Particle number (# L-1)
%      - Density (g cm-3)
%      - Velocity (m d-1)
%      - Diameter (um)
%      - Stickiness
%      - Porosity
%      - Solid volume (µm3)
%      - Fractal dimension
%      - Primary particle radius (µm)
%      - Depth (m)
%      - Organic carbon, TEP-C, opal, calcite, clay (mol)
%    For imaging depth levels, averages are categorized into:
%      - Size classifications: All, small, and large
%      - Velocity classifications: All, slow, and fast
% - Particle statistics: 
%      - Total particle count at imaging depths
%      - Largest size class, fastest velocity class
%      - Fraction of large and fast particles
% - Sources-minus-sinks (SMS) terms (mol m-3 d-1 and mol m-2 yr-1)
% - Auxiliary terms (various units)
% - Cluster and particle evolution: 
%      - Time series of the number of clusters and particles
%      - Snapshot of final clusters
%
%   Written by A. Rufas, University of Oxford
%   Version 1.0 - Completed 26 Feb 2025   
%
% Folder paths:
  % fullpathModelRunsDir      = './tests/test96/modelruns/';
  % fullpathModelInputDataDir = './tests/test96/modelinputdata/';
%
% Filenames:
  % filenameOutputSlams         = 'runsoutput.mat';
  % filenameInputRunGrid        = 'grid_run.mat';
  % filenameInputNumDepthLayers = 'waterColNumDepthLayers.txt';
%
% Choices:
    % choiceTypeGridDomain = 2; 
%       1: global
%       2: local
%       3: global with zoom in factor
    % choiceTestModelEquilibrium = false; 
%       true : the model runs for long in one location
%       false: save output for every model run year
    % choiceSaveDiskSpace = true; 
%       true : .F90 files run with #ifndef SAVE_DISK_SPACE (true when
%       running with the optimiser)
%       false: 
    % choiceSeasonalOutput = false; 
%        true : save seasonal data (as well as monthly and annual data)
%        false: only monthly and annual data
% 
% addpath(fullpathModelRunsDir);
% addpath(fullpathModelInputDataDir);
% addpath(genpath('./modelresources/'));
% addpath(genpath('./code/matlab/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Log file to record progress when running this script
logID = fopen(fullfile(fullpathModelRunsDir,'logOutputAnalysis'),'w'); % this will delete the existing file contents and create a new one

% Load configuration parameters used in the model
config = loadModelConfigurationParameters(fullpathModelInputDataDir,...
    filenameInputRunGrid,filenameInputNumDepthLayers,choiceTypeGridDomain);
    
% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - OUTPUT ARRAYS TAHT WILL BE GENERATED IN THIS SCRIPT
% -------------------------------------------------------------------------

output = struct(); % initialise the main output structure

% Sediment trap and imaging systems deployment depths
output.deployment = struct();
output.deployment.sedTrap        = NaN(config.maxNoSedTrapDeployDepths, config.nLocs);
output.deployment.imagSys        = NaN(config.maxNoImagSysDeployDepths, config.nLocs);
output.deployment.layerMidDepth  = NaN(config.maxNoDepthLayers, config.nLocs);
output.deployment.layerThickness = NaN(config.maxNoDepthLayers, config.nLocs);

% Sources-minus-sinks terms and auxiliary variables
output.sms = struct();
output.sms.monthly           = NaN(config.maxNoSMSterms, config.maxNoDepthLayers, 12, config.nLocs);
output.sms.monthlyIntegrated = NaN(config.maxNoSMSterms, 12, config.nLocs);
output.sms.annual            = NaN(config.maxNoSMSterms, config.maxNoDepthLayers, config.nLocs);
output.aux = struct();
output.aux.annual            = NaN(config.maxNoAuxTerms, config.maxNoDepthLayers, config.nLocs);
output.aux.monthly           = NaN(config.maxNoAuxTerms, config.maxNoDepthLayers, 12, config.nLocs); % memory-intensive
output.aux.seasonal          = NaN(config.maxNoAuxTerms, config.maxNoDepthLayers, 4, config.nLocs);
% if ~choiceSaveDiskSpace
% 
% end

% Water-column and seabed fluxes (1=orgC, 2=tepC, 3=calc, 4=opal, 5=clay)
output.flux = struct();
output.flux.monthly          = NaN(5, config.maxNoSedTrapDeployDepths, 12, config.nLocs); 
output.flux.annual           = NaN(5, config.maxNoSedTrapDeployDepths, 2, config.nLocs);  % 3rd dim: 1=mean, 2=std
if ~choiceSaveDiskSpace
    output.flux.seabedAnnual = zeros(3, 2, 2, config.nLocs);
end

% Particle average attributes
output.particle = struct();
output.particle.avgAttInScMonthly     = NaN(config.nDiameterClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, 12, config.nLocs);
output.particle.avgAttInVcMonthly     = NaN(config.nVelocityClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, 12, config.nLocs);
output.particle.avgAttInScAnnual      = NaN(config.nDiameterClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, config.nLocs);
output.particle.avgAttScAtDepthAnnual = NaN(3, config.maxNoImagSysDeployDepths, config.nParticleAttributes + 1, config.nLocs); % +1 for avg size
output.particle.avgAttInTyAnnual      = NaN(config.nParticleTypes, config.maxNoImagSysDeployDepths, 5, config.nLocs);
if ~choiceSaveDiskSpace
    output.particle.avgAttInVcAnnual        = NaN(config.nVelocityClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, config.nLocs);
    output.particle.avgAttVcAtDepthAnnual   = NaN(3, config.maxNoImagSysDeployDepths, config.nParticleAttributes + 1, config.nLocs); % +1 for avg velocity
    output.particle.fracMassAssocToTyAnnual = NaN(config.nParticleTypes, config.maxNoImagSysDeployDepths, config.nLocs);
    output.particle.fracPocAssocToTyAnnual  = NaN(config.nParticleTypes, config.maxNoImagSysDeployDepths, config.nLocs);
    output.particle.sizeAndVeloStatsAnnual  = NaN(5, config.maxNoSedTrapDeployDepths, config.nLocs); % 1=particle number, 2=largest size class, 3=fastest velo class, 4=frac large particles, 5=frac fast particles
end
if ~choiceSaveDiskSpace 
	output.particle.snapshots = []; % for out_cluster_i.bin, out_cluster_r.bin
end

output.seasonal.particle.avgAttInSc        = NaN(config.nDiameterClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, 4, config.nLocs);
output.seasonal.particle.avgAttInVc        = NaN(config.nVelocityClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, 4, config.nLocs);
output.seasonal.particle.avgAttScAtDepth   = NaN(3, config.maxNoImagSysDeployDepths, config.nParticleAttributes + 1, 4, config.nLocs); % +1 for avg size
output.seasonal.particle.avgAttVcAtDepth   = NaN(3, config.maxNoImagSysDeployDepths, config.nParticleAttributes + 1, 4, config.nLocs); % +1 for avg velocity

% Seasonal calculations
if choiceSeasonalOutput
    output.seasonal = struct();
    
    % Sources-minus-sinks
    output.seasonal.sms = NaN(config.maxNoSMSterms, config.maxNoDepthLayers, 4, config.nLocs);
    
    % Fluxes
    output.seasonal.flux       = NaN(5, config.maxNoSedTrapDeployDepths, 4, config.nLocs); % 1=orgC, 2=tepC, 3=calc, 4=opal, 5=clay
    output.seasonal.fluxSeabed = NaN(3, 4, 2, config.nLocs); % 1=orgC+tepC, 2=calc, 3=opal / 1=absolute value, 2=fraction
    
    % Particle attributes
    output.seasonal.particle.avgAttInSc        = NaN(config.nDiameterClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, 4, config.nLocs);
    output.seasonal.particle.avgAttInVc        = NaN(config.nVelocityClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes, 4, config.nLocs);
    output.seasonal.particle.avgAttScAtDepth   = NaN(3, config.maxNoImagSysDeployDepths, config.nParticleAttributes + 1, 4, config.nLocs); % +1 for avg size
    output.seasonal.particle.avgAttVcAtDepth   = NaN(3, config.maxNoImagSysDeployDepths, config.nParticleAttributes + 1, 4, config.nLocs); % +1 for avg velocity
    output.seasonal.particle.avgAttInTy        = NaN(config.nParticleTypes, config.maxNoImagSysDeployDepths, 5, 4, config.nLocs);
    output.seasonal.particle.fracMassAssocToTy = NaN(config.nParticleTypes, config.maxNoImagSysDeployDepths, 4, config.nLocs);
    output.seasonal.particle.fracPocAssocToTy  = NaN(config.nParticleTypes, config.maxNoImagSysDeployDepths, 4, config.nLocs);
    output.seasonal.particle.sizeAndVeloStats  = NaN(5, config.maxNoSedTrapDeployDepths, 4, config.nLocs); % 1=particle number, 2=largest size class, 3=fastest velo class, 4=frac large particles, 5=frac fast particles
end
       
% If we want to save output for every model run year
if choiceTestModelEquilibrium
    maxNoTimeStepsRecording = 12*config.nYears; % nMonths * N years (change the number of years accordingly)
    output.evolution = struct();
    output.evolution.amountOfParticles          = NaN(maxNoTimeStepsRecording, 2, config.nLocs); % 1=num of clusters, 2=num of particles
    output.evolution.flux                       = NaN(5, config.maxNoSedTrapDeployDepths, maxNoTimeStepsRecording, config.nLocs);
    output.evolution.fluxSeasonal               = NaN(5, config.maxNoSedTrapDeployDepths, 4, config.nYears, config.nLocs);
    output.evolution.fluxAnnual                 = NaN(5, config.maxNoSedTrapDeployDepths, maxNoTimeStepsRecording, config.noYears, config.nLocs);
    output.evolution.particleAvgAttInScSeasonal = NaN(config.nDiameterClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes-1, 4, config.nYears, config.nLocs);    
    output.evolution.particleAvgAttInScAnnual   = NaN(config.nDiameterClasses, config.maxNoImagSysDeployDepths, config.nParticleAttributes-1, config.nYears, config.nLocs);
    output.evolution.particleNumberInScAnnual   = NaN(config.nDiameterClasses, config.maxNoImagSysDeployDepths, config.nYears, config.nLocs);
end

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 7 - READ SLAMS OUTPUT FILES
% -------------------------------------------------------------------------

tic; % start timing the processing

for iRun = 1:config.nLocs

    fprintf(logID, 'Location id: %d\n', iRun);
    fullpathThisRunDir = fullfile(fullpathModelRunsDir, sprintf('run_%d', iRun));

    % Define the full paths to the files that will be read by this script
    runStatus                     = checkStatus(fullpathThisRunDir,config.filenameOutputStatus,logID,iRun);
    controlFile                   = fullfile(fullpathThisRunDir,config.filenameOutputControl);
    zubFile                       = fullfile(fullpathThisRunDir,config.filenameModelZub);
    zlbFile                       = fullfile(fullpathThisRunDir,config.filenameModelZlb);
    evolNumClustersFile           = fullfile(fullpathThisRunDir,config.filenameOutputNumClusters);
    evolNumParticlesFile          = fullfile(fullpathThisRunDir,config.filenameOutputNumParticles);
    clusterIntAttributesFile     = fullfile(fullpathThisRunDir,config.filenameOutputClusterInteger);
    clusterRealAttributesFile     = fullfile(fullpathThisRunDir,config.filenameOutputClusterReal);
    smsFile                       = fullfile(fullpathThisRunDir,config.filenameOutputSms);
    auxFile                       = fullfile(fullpathThisRunDir,config.filenameOutputAux);
    smsIntegratedFile             = fullfile(fullpathThisRunDir,config.filenameOutputSmsIntegrated);
    fluxFile                      = fullfile(fullpathThisRunDir,config.filenameOutputFlux);
    fluxSfFile                    = fullfile(fullpathThisRunDir,config.filenameOutputFluxSf);
    particleAttsScFile            = fullfile(fullpathThisRunDir,config.filenameOutputAvgAttsSizeClass);
    particleAttsScSfFile          = fullfile(fullpathThisRunDir,config.filenameOutputAvgAttsSizeClassSf);
    particleAttsTyFile            = fullfile(fullpathThisRunDir,config.filenameOutputAvgAttsParticleType);
    particleAttsTySfFile          = fullfile(fullpathThisRunDir,config.filenameOutputAvgAttsParticleTypeSf);
    particleAttsVcFile            = fullfile(fullpathThisRunDir,config.filenameOutputAvgAttsVeloClass);
    particleAttsVcSfFile          = fullfile(fullpathThisRunDir,config.filenameOutputAvgAttsVeloClassSf);

    if strcmp(runStatus, 'completed')
        
        %------------------------------------------------------------------
        % GENERAL INFO
        %------------------------------------------------------------------

        % Read control file
        if isfile(controlFile)
            control = FNread_binary(controlFile,[9,1],'*int32','ieee-be');
            nTimeSteps                = control(1);
            nLocalSedTrapDeployDepths = control(4);
            nLocalImagSysDeployDepths = control(5);
            nWriteStepsSedTrap        = control(3); % SedTrap, AvgAtt and Stats are the same
            nWriteStepsAvgAtt         = control(6);
            nWriteStepsStats          = control(7);
            nWriteStepsEulerian       = nWriteStepsSedTrap;
            nWriteStepsSnapshots      = control(9);
            %NUM_YEARS = nTimeSteps/NUM_TIMESTEPS_PER_YEAR;
        else
            fprintf(logID, 'Control file "%s" does not exist.', controlFile);
            return; % exit the script or function
        end

        % .................................................................
         
        % Read depths
        nLocalDepthLayersOriginal = sum(~isnan(config.gridDepths(:,iRun))); % this is from the grid array Zsb
        nLocalDepthLayersActual = config.nDepthLayersUsed(iRun); % this is from 'waterColNumDepthLayers.txt', which could have been modified before running locally

        if isfile(zubFile) && isfile(zlbFile)
            layerUpperBoundary = FNread_binary(zubFile,[nLocalDepthLayersOriginal, 1],'float64','ieee-be'); 
            layerLowerBoundary = FNread_binary(zlbFile,[nLocalDepthLayersOriginal, 1],'float64','ieee-be');
        else
           fprintf(logID, 'Depth files "%s" or "%s" do not exist.', zubFile, zlbFile);
           return; % exit the script or function 
        end

        if (nLocalDepthLayersActual ~= nLocalDepthLayersOriginal)
            layerUpperBoundary = layerUpperBoundary(1:nLocalDepthLayersActual);
            layerLowerBoundary = layerLowerBoundary(1:nLocalDepthLayersActual);
        end
      
        layerThickness = layerLowerBoundary - layerUpperBoundary; % m
        layerMidDepth = layerLowerBoundary - layerThickness * 0.5;
        seafloorDepth = layerLowerBoundary(nLocalDepthLayersActual);

        % Depths at which the model is collecting particles with sediment traps
        sedTrapWaterColDeployDepthsLocal = zeros(nLocalSedTrapDeployDepths,1);
        for i = 1:nLocalSedTrapDeployDepths
            sedTrapWaterColDeployDepthsLocal(i) = config.availSedTrapDeployDepths(i);
        end
        sedTrapAllDeployDepthsLocal = cat(1, sedTrapWaterColDeployDepthsLocal(:), seafloorDepth);
        nLastDeployDepthSedTrapLocal = numel(sedTrapAllDeployDepthsLocal);
        
		% Depths at which the model is imaging particles
        imagSysWaterColDeployDepthsLocal = zeros(nLocalImagSysDeployDepths,1);
        for i = 1:nLocalImagSysDeployDepths
            imagSysWaterColDeployDepthsLocal(i) = config.availImagSysDeployDepths(i);
        end
        imagSysAllDeployDepthsLocal = cat(1,imagSysWaterColDeployDepthsLocal(:),seafloorDepth);
        nLastDeployDepthImagSysLocal = numel(imagSysAllDeployDepthsLocal);
        
        % Store depth-related data in the output structure
        output.deployment.layerThickness(1:nLocalDepthLayersActual,iRun) = layerThickness;
        output.deployment.layerMidDepth(1:nLocalDepthLayersActual,iRun) = layerMidDepth;
        output.deployment.sedTrap(1:nLastDeployDepthSedTrapLocal,iRun) = sedTrapAllDeployDepthsLocal(:);
        output.deployment.imagSys(1:nLastDeployDepthImagSysLocal,iRun) = imagSysAllDeployDepthsLocal(:);

        %------------------------------------------------------------------
        % LOAD PARTICLE ARRAY DATA
        %------------------------------------------------------------------
        
        if ~choiceSaveDiskSpace
            
            if isfile(evolNumClustersFile) && isfile(evolNumParticlesFile)
                % Read no. clusters and particles at every time step
                nClustersEvol = FNread_binary(evolNumClustersFile,[nWriteStepsSedTrap, 1],'*int32','ieee-be'); 
                nParticlesEvol = FNread_binary(evolNumParticlesFile,[nWriteStepsSedTrap, 1],'*int64','ieee-be'); 
            
                % Track the last values for reporting
                nClusters = nClustersEvol(end);
                nParticles = nParticlesEvol(end);
                
                if choiceTestModelEquilibrium
                    output.evolution.amountOfParticles(1:nWriteStepsSedTrap,1,iRun) = nClustersEvol;
                    output.evolution.amountOfParticles(1:nWriteStepsSedTrap,2,iRun) = nParticlesEvol;
                end
            else
                fprintf(logID, 'Error: Particle evolution files "%s" or "%s" not found.\n', ...
                    evolNumClustersFile, evolNumParticlesFile);
                return;
            end
    
        	if isfile(clusterIntAttributesFile) && isfile(clusterRealAttributesFile)
            	clusterInt = FNread_binary(clusterIntAttributesFile, [config.maxNoClusters, 6, nWriteStepsSnapshots],'*int32','ieee-be'); % integer attributes
            	clusterReal = FNread_binary(clusterRealAttributesFile, [config.maxNoClusters, 19, nWriteStepsSnapshots],'float64','ieee-be'); % real attributes
            	output.particle.snapshots = NaN(config.maxNoClustersCutoff,26,nWriteStepsSnapshots,config.nLocs);

            	% Store cluster data in the structure
            	output.particle.snapshots(:,1:6,:,iRun) = clusterInt(1:config.maxNoClustersCutoff,:,:);
            	output.particle.snapshots(:,7:25,:,iRun) = clusterReal(1:config.maxNoClustersCutoff,:,:);
        	else
            	fprintf(logID, 'Error: Cluster files "%s" or "%s" not found.\n', ...
                	clusterIntAttributesFile, clusterRealAttributesFile);
            	return;
        	end
        
        end
        
        %------------------------------------------------------------------
        % LOAD SOURCES-MINUS-SINKS AND AUXILLIARY DATA DATA
        %------------------------------------------------------------------

        if isfile(smsFile) && isfile(auxFile)
            modelSms = FNread_binary(smsFile, [config.maxNoSMSterms,...
                nLocalDepthLayersActual, nWriteStepsEulerian],'float64','ieee-be');
            modelAux = FNread_binary(auxFile, [config.maxNoAuxTerms,...
                nLocalDepthLayersActual, nWriteStepsEulerian],'float64','ieee-be'); 

            if isfile(smsIntegratedFile)
                if ~choiceSaveDiskSpace
                    modelSmsIntegrated = FNread_binary(smsIntegratedFile,...
                        [config.maxNoSMSterms, nWriteStepsEulerian],'float64','ieee-be');
                else
                    modelSmsIntegrated = [];
                end
            end
        else
           fprintf(logID, 'SMS and aux. files "%s" or "%s" do not exist.', smsFile, auxFile);
           return; 
        end
        
        %------------------------------------------------------------------
        % LOAD FLUX DATA
        %------------------------------------------------------------------
        
        % 1=orgC, 2=TEP-C, 3=calcite, 4=opal, 5=clay, mg m-2 d-1
        
        if isfile(fluxFile) && isfile(fluxSfFile)
            modelFluxWaterCol = FNread_binary(fluxFile, [5,...
                nLocalSedTrapDeployDepths, nWriteStepsSedTrap],'float64','ieee-be'); 
            modelFluxSeafloor = FNread_binary(fluxSfFile, [5, 1,...
                nWriteStepsSedTrap],'float64','ieee-be'); 

            % Apply flux detection limit
            modelFluxWaterCol(modelFluxWaterCol < config.fluxDetectionLimit) = 0;
            modelFluxSeafloor(modelFluxSeafloor < config.fluxDetectionLimit) = 0;

            % Combine water column and seafloor flux data
            modelFlux = cat(2,modelFluxWaterCol,modelFluxSeafloor);

            % If testing model equilibrium, store the evolution of flux
            if choiceTestModelEquilibrium
                output.evolution.flux(:,1:(nLocalSedTrapDeployDepths + 1),1:nWriteStepsSedTrap,iRun) = modelFlux;
            end
        else
            fprintf(logID, 'Error: Flux data files "%s" or "%s" not found.\n', fluxFile, fluxSfFile);
            return;
        end

        %------------------------------------------------------------------
        % LOAD PARTICLE SPECTRA DATA
        %------------------------------------------------------------------

        % For the water column particles, the units are # particles L-1,
        % and for the seafloor particles, the units are # particles.

        if isfile(particleAttsScFile) && isfile(particleAttsScSfFile)
            modelAvgAttSizeSpec = FNread_binary(particleAttsScFile,...
                [config.nDiameterClasses, nLocalImagSysDeployDepths,...
                config.nParticleAttributes, nWriteStepsAvgAtt], 'float64', 'ieee-be');
            modelAvgAttSizeSpecSf = FNread_binary(particleAttsScSfFile,...
                [config.nDiameterClasses, 1, config.nParticleAttributes,...
                nWriteStepsAvgAtt], 'float64', 'ieee-be');
        else
            fprintf(logID, 'Error: particle avg. attributes data files "%s" or "%s" not found.\n',...
               particleAttsScFile, particleAttsScFile);
            return;
        end
        if isfile(particleAttsTyFile) && isfile(particleAttsTySfFile)
            modelAvgAttParticleType = FNread_binary(particleAttsTyFile,...
                [config.nParticleTypes, nLocalImagSysDeployDepths, 7, nWriteStepsAvgAtt], 'float64', 'ieee-be');
            modelAvgAttParticleTypeSf = FNread_binary(particleAttsTySfFile,...
                [config.nParticleTypes, 1, 7, nWriteStepsAvgAtt], 'float64', 'ieee-be');
        else
            fprintf(logID, 'Error: particle avg. attributes data files "%s" or "%s" not found.\n',...
               particleAttsTyFile, particleAttsTySfFile);
            return;
        end
        if isfile(particleAttsVcFile) && isfile(particleAttsVcSfFile)
            modelAvgAttVeloSpec = FNread_binary(particleAttsVcFile,...
                [config.nVelocityClasses, nLocalImagSysDeployDepths,...
                config.nParticleAttributes, nWriteStepsAvgAtt], 'float64', 'ieee-be');
            modelAvgAttVeloSpecSf = FNread_binary(particleAttsVcSfFile,...
                [config.nVelocityClasses, 1, config.nParticleAttributes,...
                nWriteStepsAvgAtt], 'float64', 'ieee-be');
        else
            fprintf(logID, 'Error: particle avg. attributes data files "%s" or "%s" not found.\n',...
               particleAttsVcFile, particleAttsVcSfFile);
            return;
        end
        
        %-------------------------------------------------------------------
        % PROCESS FLUXES AND PARTICLE DATA BY SIZE CLASS FOR ALL MODEL RUN YEARS
		%-------------------------------------------------------------------    
      
        if choiceTestModelEquilibrium
            output = calculateAverageDataForAllRunYears(output,iRun,config,...
                modelFlux,modelAvgAttSizeSpec,modelAvgAttSizeSpecSf,...
                nLastDeployDepthSedTrapLocal,nLastDeployDepthImagSysLocal,...
                nLocalImagSysDeployDepths); 
        end
       
        % The following sections consider only the last model run year,
        % which is considered to be "at equilibrium"
        
        %------------------------------------------------------------------
        % PROCESS SOURCES-MINUS-SINKS AND AUXILLIARY TERMS AT LAST MODEL RUN YEAR
        %------------------------------------------------------------------

        output = processSourcesMinusSinksAndAuxilliaryTerms(output,iRun,config,...
            modelSms,modelSmsIntegrated,modelAux,nWriteStepsEulerian,...
            nLocalDepthLayersActual,choiceSeasonalOutput,choiceSaveDiskSpace);

        %------------------------------------------------------------------
        % PROCESS FLUXES AND PARTICLE DATA BY SIZE CLASS, VELOCITY CLASS AND
        % TYPE AT LAST MODEL RUN YEAR 
		%------------------------------------------------------------------

        output = processEquilibriumFluxAndParticleAttributes(output,iRun,config,...
            modelFlux,modelAvgAttSizeSpec,modelAvgAttSizeSpecSf,modelAvgAttParticleType,...
            modelAvgAttParticleTypeSf,modelAvgAttVeloSpec,modelAvgAttVeloSpecSf,...
            nLastDeployDepthSedTrapLocal,nLastDeployDepthImagSysLocal,nLocalImagSysDeployDepths,...
            nWriteStepsSedTrap,nWriteStepsAvgAtt,choiceSeasonalOutput,choiceSaveDiskSpace);

        %------------------------------------------------------------------
        % PROCESS SEAFLOOR MAP OF MATERIALS AT LAST MODEL RUN YEAR
        %------------------------------------------------------------------

        output = processSeafloorFluxAverages(output,iRun,config,nLastDeployDepthSedTrapLocal,...
           choiceSeasonalOutput,choiceSaveDiskSpace);

        %------------------------------------------------------------------
		% PROCESS AVERAGE PARTICLE DATA BY DEPTH AND LARGE/SMALL & SLOW/FAST
		% CLASS AT LAST MODEL RUN YEAR
        %------------------------------------------------------------------

        output = processParticleDataBySizeAndVelocityClass(output,iRun,config,...
            nLastDeployDepthImagSysLocal,choiceSeasonalOutput,choiceSaveDiskSpace);

    end % end checking status of run
end % iRun

% Finalise and close files
vars = who('output*');
save(fullfile(fullpathModelRunsDir,filenameOutputSlams),vars{:},'-v7.3'); % this will delete the existing file contents and create a new one
fprintf(logID, 'Processing completed in %.2f seconds.\n', toc);    
fclose(logID);

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

% *************************************************************************

function status = checkStatus(fullpathThisRunDir,filenameOutputStatus,logID,iRun)

    % Construct the full path to the status file
    statusFile = fullfile(fullpathThisRunDir,filenameOutputStatus);
    
    % Check if the file exists
    if ~isfile(statusFile)
        throwError(logID,'The file "%s" does not exist in the specified folder %d',filenameOutputStatus,iRun);
    end

    % Open and read the status file
    fileID = fopen(statusFile, 'r');
    status = strtrim(fgetl(fileID)); % read the first line and remove extra whitespace
    fclose(fileID);

    % Display message based on the status
    switch lower(status)
        case 'completed'
            fprintf(logID, 'Status is "completed". Proceeding with array reading...');
            
        case 'no clusters seeded'
            fprintf(logID, 'Status is "no clusters seeded". No further action will be taken.');
            
        case 'water too shallow'
            fprintf(logID, 'Status is "water too shallow". No further action will be taken.');
            
        case 'in progress'
            fprintf(logID, 'Status is "in progress". Waiting for completion.');
            
        otherwise
            throwError(logID, 'Unexpected status found in "%s": %s', filenameOutputStatus, status);
    end

end % checkStatus

% *************************************************************************

function weightedAverage = computeWeightedAverage(data,weights,sumDim)

    % Add the values of the months in the same season and average them 
    % to get seasonal aux. To do this properly, we need to calculate 
    % the weighted average considering that each month has a different 
    % number of days
        
    % Sum the weighted values along the specified axis
    weightedSum = sum(data .* reshape(weights, 1, 1, []), sumDim);
    
    % Divide by the total weights
    totalWeight = sum(weights);
    weightedAverage = weightedSum / totalWeight;
    
end % weightedAverage

% *************************************************************************

function seasonalData = computeSeasonalData(monthlyData,config)

    % Helper function for weighted seasonal aggregation of data.
    
    % Add the values of the months in the same season and average them 
    % to get seasonal values. To do this properly, we need to calculate 
    % the weighted average considering that each month has a different 
    % number of days. 

    nDims = ndims(monthlyData); % determine if 'monthlyData' includes the year dimension
    
    % Handle the case where monthlyData has the year dimension or not
    if nDims == 4
        nYears = size(monthlyData,4); % get the number of years
        seasonalData = zeros(size(monthlyData,1),size(monthlyData,2),nSeasons,nYears);
        
        for iSeason = 1:config.nSeasons
            monthsInSeason = (3 * (iSeason - 1) + 1):(3 * iSeason); % define months for the season
            seasonalData(:,:,iSeason,:) = sum(...
                monthlyData(:,:,monthsInSeason,:) .* reshape(config.nDaysPerMonth(monthsInSeason),1,1,[],1), 3)...
                ./ config.nDaysPerSeason(iSeason);
        end
        
    elseif nDims == 3
        seasonalData = zeros(size(monthlyData,1),size(monthlyData,2),config.nSeasons);
        
        for iSeason = 1:config.nSeasons
            monthsInSeason = (3 * (iSeason - 1) + 1):(3 * iSeason); % define months for the season
            seasonalData(:,:,iSeason) = sum(...
                monthlyData(:,:,monthsInSeason) .* reshape(config.nDaysPerMonth(monthsInSeason),1,1,[]), 3)...
                ./ config.nDaysPerSeason(iSeason);
        end

    end
    
end % computeSeasonalData

% *************************************************************************

function annualData = computeAnnualData(monthlyData,config)

    % Helper function for weighted annual aggregation of data.
    
    % Compute total number of days in the year
    totalDays = sum(config.nDaysPerMonth);
    
    % Reshape nDaysPerMonth to match the dimensionality of monthlyData
    nDims = ndims(monthlyData);
    
    if nDims == 4
        % If monthlyData is 4D, reshape nDaysPerMonth as (1, 1, 12, 1)
        reshapedDays = reshape(config.nDaysPerMonth, 1, 1, [], 1);
    elseif nDims == 3
        % If monthlyData is 3D, reshape nDaysPerMonth as (1, 1, 12)
        reshapedDays = reshape(config.nDaysPerMonth, 1, 1, []);
    end
    
    % Calculate annual data by summing the weighted monthly data and 
    % normalising by total days
    annualData = sum(monthlyData .* reshapedDays, 3) ./ totalDays;
    
end % computeAnnualData

% *************************************************************************

function classRange = getClassRange(classIdx,config)

    % Helper function to get class ranges based on thresholds

    switch classIdx
        case 1 % All Sc
            classRange = 1:config.nDiameterClasses;
        case 2 % Small Sc
            classRange = 1:config.idxLargeScThreshold;
        case 3 % Large Sc
            classRange = (config.idxLargeScThreshold + 1):config.nDiameterClasses;
        case 4 % All Vc
            classRange = 1:config.nVelocityClasses;
        case 5 % Slow Vc
            classRange = 1:config.idxFastVcThreshold;
        case 6 % Fast Vc
            classRange = (config.idxFastVcThreshold + 1):config.nVelocityClasses;
    end
    
end % getClassRange

% *************************************************************************

function outputArray = aggregateParticleData(timeMode,monthlyData,weightField,...
    config)
        
    % For the particle attributes, average as in the output files from SLAMS,
    % that is, multiply each attribute value by the number of particles
    % it represents, sum all (att x pNum) and divide by the total number of
    % particles. This way, we get a pondered attribute value. For the
    % particle number, just average when there are particle instances.
    
    % In order to average the number of particles, we divide the monthly 
    % pondered adddition by the sum of the number of days in the season, 
    % e.g. (1000*31 + 800*28 + 950*30)/(31+28+30)
    
    % In order to average the attributes, we divide the pondered addition
    % by the number of particles (i.e., 2.8*1000*31 + 2.7*800*28 + 2.5*950*30)/(1000*31+800*28+950*30))
              
    switch lower(timeMode)
        case 'annual'
            nPeriods = 1;
            periodMapping = ones(1, 12); % all months belong to the same annual period
            totalDaysPerPeriod = sum(config.nDaysPerMonth); % total days in the year
        case 'seasonal'
            nPeriods = 4;
            periodMapping = ceil((1:12) / 3); % map months to 4 seasons
            totalDaysPerPeriod = config.nDaysPerSeason;
    end

    dataArraySize = size(monthlyData);
    weightsArraySize = size(weightField);
    
    % Adjust output size based on whether it's 3D or 4D
    if length(dataArraySize) == 3 % particle data without attribute dimension (e.g., pNum)
        outputArray = zeros(dataArraySize(1), dataArraySize(2), nPeriods);
    elseif length(dataArraySize) == 4 % particle data with attribute dimension (e.g., avgAtt)
        outputArray = zeros(dataArraySize(1), dataArraySize(2), dataArraySize(3), nPeriods);
        weightFieldAdd = zeros(weightsArraySize(1), weightsArraySize(2), nPeriods);
    end
    
    % Aggregate data
    for iMonth = 1:12
        iPeriod = periodMapping(iMonth);   
        if length(dataArraySize) == 3
            outputArray(:,:,iPeriod) = outputArray(:,:,iPeriod) + ...
                monthlyData(:,:,iMonth) .* config.nDaysPerMonth(iMonth);
        elseif length(dataArraySize) == 4
            for iAtt = 1:dataArraySize(3)
                outputArray(:,:,iAtt,iPeriod) = outputArray(:,:,iAtt,iPeriod) + ...
                    monthlyData(:,:,iAtt,iMonth) .* weightField(:,:,iMonth) .* config.nDaysPerMonth(iMonth);
                weightFieldAdd(:,:,iPeriod) = weightFieldAdd(:,:,iPeriod) + ...
                    weightField(:,:,iMonth) .* config.nDaysPerMonth(iMonth);
            end
        end
    end

    % Finalise aggregation by dividing sums
    for iPeriod = 1:nPeriods
        if length(dataArraySize) == 3
            outputArray(:,:,iPeriod) = outputArray(:,:,iPeriod) / totalDaysPerPeriod(iPeriod);
        elseif length(dataArraySize) == 4
            for iAtt = 1:dataArraySize(3)
                outputArray(:,:,iAtt,iPeriod) = outputArray(:,:,iAtt,iPeriod) ./ weightFieldAdd(:,:,iPeriod);
            end
        end
    end
    
end % aggregateParticleData

% *************************************************************************

function output = processEquilibriumFluxAndParticleAttributes(output,iRun,config,...
    modelFlux,modelAvgAttSizeSpec,modelAvgAttSizeSpecSf,modelAvgAttParticleType,...
    modelAvgAttParticleTypeSf,modelAvgAttVeloSpec,modelAvgAttVeloSpecSf,...
    nLastDeployDepthSedTrapLocal,nLastDeployDepthImagSysLocal,nLocalImagSysDeployDepths,...
    nWriteStepsSedTrap,nWriteStepsAvgAtt,choiceSeasonalOutput,choiceSaveDiskSpace)

    % Initialise output arrays for monthly data
    fluxMonthly              = zeros(5,nLastDeployDepthSedTrapLocal,12);
    pNumInScMonthly          = zeros(config.nDiameterClasses,nLastDeployDepthImagSysLocal,12);
    avgAttInScMonthly        = zeros(config.nDiameterClasses,nLastDeployDepthImagSysLocal,config.nParticleAttributes,12);
    pNumInTyMonthly          = zeros(config.nParticleTypes,nLastDeployDepthImagSysLocal,12);
    avgAttInTyMonthly        = zeros(config.nParticleTypes,nLastDeployDepthImagSysLocal,4,12);
    pNumInVcMonthly          = zeros(config.nVelocityClasses,nLastDeployDepthImagSysLocal,12); 
    avgAttInVcMonthly        = zeros(config.nVelocityClasses,nLastDeployDepthImagSysLocal,config.nParticleAttributes,12);
    fracMassAssocToTyMonthly = zeros(config.nParticleTypes,nLastDeployDepthImagSysLocal,12);
    fracPocAssocToTyMonthly  = zeros(config.nParticleTypes,nLastDeployDepthImagSysLocal,12);

    % .................................................................

    % Extract monthly data for the last year
    for iMonth = 1:12

        iSliceFlux = (nWriteStepsSedTrap - 12) + iMonth;
        fluxMonthly(:,:,iMonth) = modelFlux(:,:,iSliceFlux);

        iSliceAtts = (nWriteStepsAvgAtt-12) + iMonth;
        pNumInScMonthly(:,(1:nLocalImagSysDeployDepths),iMonth) = modelAvgAttSizeSpec(:,:,1,iSliceAtts);
        pNumInScMonthly(:,nLastDeployDepthImagSysLocal,iMonth) = modelAvgAttSizeSpecSf(:,1,1,iSliceAtts);
        pNumInTyMonthly(:,(1:nLocalImagSysDeployDepths),iMonth) = modelAvgAttParticleType(:,:,1,iSliceAtts);
        pNumInTyMonthly(:,nLastDeployDepthImagSysLocal,iMonth) = modelAvgAttParticleTypeSf(:,1,1,iSliceAtts);
        pNumInVcMonthly(:,(1:nLocalImagSysDeployDepths),iMonth) = modelAvgAttVeloSpec(:,:,1,iSliceAtts);
        pNumInVcMonthly(:,nLastDeployDepthImagSysLocal,iMonth) = modelAvgAttVeloSpecSf(:,1,1,iSliceAtts);

        avgAttInScMonthly(:,(1:nLocalImagSysDeployDepths),:,iMonth) = modelAvgAttSizeSpec(:,:,:,iSliceAtts);
        avgAttInScMonthly(:,nLastDeployDepthImagSysLocal,:,iMonth) = modelAvgAttSizeSpecSf(:,1,:,iSliceAtts);	
        avgAttInTyMonthly(:,(1:nLocalImagSysDeployDepths),:,iMonth) = modelAvgAttParticleType(:,:,(4:7),iSliceAtts); % avg. diameter, density, velocity, porosity
        avgAttInTyMonthly(:,nLastDeployDepthImagSysLocal,:,iMonth) = modelAvgAttParticleTypeSf(:,1,(4:7),iSliceAtts);
        avgAttInVcMonthly(:,(1:nLocalImagSysDeployDepths),:,iMonth) = modelAvgAttVeloSpec(:,:,:,iSliceAtts);
        avgAttInVcMonthly(:,nLastDeployDepthImagSysLocal,:,iMonth) = modelAvgAttVeloSpecSf(:,1,:,iSliceAtts);

        for iImagDepth = 1:nLocalImagSysDeployDepths
            totMassAtId = squeeze(sum(modelAvgAttParticleType(:,iImagDepth,2,iSliceAtts)));
            totPocAtId = squeeze(sum(modelAvgAttParticleType(:,iImagDepth,3,iSliceAtts)));
            if totMassAtId > 0
                for iTy = 1:config.nParticleTypes
                    fracMassAssocToTyMonthly(iTy,iImagDepth,iMonth) = modelAvgAttParticleType(iTy,iImagDepth,2,iSliceAtts)./totMassAtId;
                    fracPocAssocToTyMonthly(iTy,iImagDepth,iMonth) = modelAvgAttParticleType(iTy,iImagDepth,3,iSliceAtts)./totPocAtId;
                end
            end
        end

        totMassAtSf = squeeze(sum(modelAvgAttParticleTypeSf(:,1,2,iSliceAtts)));
        totPocAtSf = squeeze(sum(modelAvgAttParticleTypeSf(:,1,3,iSliceAtts)));
        if totMassAtSf > 0
            for iTy = 1:config.nParticleTypes
                fracMassAssocToTyMonthly(:,nLastDeployDepthImagSysLocal,iMonth) = modelAvgAttParticleTypeSf(iTy,1,2,iSliceAtts)./totMassAtSf;
                if totPocAtSf > 0
                    fracPocAssocToTyMonthly(:,nLastDeployDepthImagSysLocal,iMonth) = modelAvgAttParticleTypeSf(iTy,1,3,iSliceAtts)./totPocAtSf;
                end
            end
        end

    end % iMonth 

    % Save to output structure 
    output.flux.monthly(:,1:nLastDeployDepthSedTrapLocal,:,iRun) = fluxMonthly;
    output.particle.avgAttInScMonthly(:,1:nLastDeployDepthImagSysLocal,1,:,iRun) = pNumInScMonthly; % particle number per Sc 
    output.particle.avgAttInScMonthly(:,1:nLastDeployDepthImagSysLocal,2:config.nParticleAttributes,:,iRun) =...
        avgAttInScMonthly(:,:,(2:config.nParticleAttributes),:); % attributes (the 1st attribute is discarded as it's the pNum transformed) 
    output.particle.avgAttInVcMonthly(:,1:nLastDeployDepthImagSysLocal,1,:,iRun) = pNumInVcMonthly; % particle number per Vc
    output.particle.avgAttInVcMonthly(:,1:nLastDeployDepthImagSysLocal,2:config.nParticleAttributes,:,iRun) =...    
        avgAttInVcMonthly(:,:,(2:config.nParticleAttributes),:); % attributes (the 1st attribute is discarded as it's the pNum transformed) 
    
    % .................................................................

    % Calculate seasonal and annual values for flux

    if choiceSeasonalOutput
        fluxSeasonal = computeSeasonalData(fluxMonthly,config);
    end
    fluxAnnual = computeAnnualData(fluxMonthly,config);

    % Save to output structure
    if choiceSeasonalOutput
        output.seasonal.flux(:,1:nLastDeployDepthSedTrapLocal,:,iRun) = fluxSeasonal;
    end
    output.flux.annual(:,1:nLastDeployDepthSedTrapLocal,1,iRun) = fluxAnnual(:,:); 
    output.flux.annual(:,1:nLastDeployDepthSedTrapLocal,2,iRun) = std(fluxMonthly(:,:,:),config.nDaysPerMonth(:),3);

    % .................................................................

    % Calculate seasonal and annual values for particle numbers

    if choiceSeasonalOutput
        pNumInScSeasonal = computeSeasonalData(pNumInScMonthly,config);
        pNumInTySeasonal = computeSeasonalData(pNumInTyMonthly,config);
        if ~choiceSaveDiskSpace
            pNumInVcSeasonal = computeSeasonalData(pNumInVcMonthly,config);
        end

        % Handle NaN
        pNumInScSeasonal(isnan(pNumInScSeasonal)) = 0;
        pNumInTySeasonal(isnan(pNumInTySeasonal)) = 0;
        if ~choiceSaveDiskSpace
            pNumInVcSeasonal(isnan(pNumInVcSeasonal)) = 0;
        end
    end

    pNumInScAnnual = computeAnnualData(pNumInScMonthly,config);
    pNumInTyAnnual = computeAnnualData(pNumInTyMonthly,config);
    if ~choiceSaveDiskSpace
        pNumInVcAnnual = computeAnnualData(pNumInVcMonthly,config);
    end

    % Handle NaN
    pNumInScAnnual(isnan(pNumInScAnnual)) = 0;
    pNumInTyAnnual(isnan(pNumInTyAnnual)) = 0;
    if ~choiceSaveDiskSpace
        pNumInVcAnnual(isnan(pNumInVcAnnual)) = 0;
    end

    % .................................................................

    % Calculate seasonal and annual values for particle attributes

    if choiceSeasonalOutput
        % Process seasonal
        avgAttInScSeasonal = aggregateParticleData('seasonal',avgAttInScMonthly,pNumInScMonthly,config);
        avgAttInTySeasonal = aggregateParticleData('seasonal',avgAttInTyMonthly,pNumInTyMonthly,config);
        if ~choiceSaveDiskSpace
            avgAttInVcSeasonal = aggregateParticleData('seasonal',avgAttInVcMonthly,...
                pNumInVcMonthly,config);
        end
        fracMassAssocToTySeasonal = aggregateParticleData('seasonal',fracMassAssocToTyMonthly,...
            pNumInTyMonthly,config);
        fracPocAssocToTySeasonal = aggregateParticleData('seasonal',fracPocAssocToTyMonthly,...
            pNumInTyMonthly,config);

        % Handle NaN
        avgAttInScSeasonal(isnan(avgAttInScSeasonal)) = 0;
        avgAttInTySeasonal(isnan(avgAttInTySeasonal)) = 0;
        if ~choiceSaveDiskSpace
            avgAttInVcSeasonal(isnan(avgAttInVcSeasonal)) = 0;
        end
        fracMassAssocToTySeasonal(isnan(fracMassAssocToTySeasonal)) = 0;
        fracPocAssocToTySeasonal(isnan(fracPocAssocToTySeasonal)) = 0;

        % Save to output structure
        output.seasonal.particle.avgAttInSc(:,1:nLastDeployDepthImagSysLocal,1,:,iRun) =...
            pNumInScSeasonal; % particle number per Sc
        output.seasonal.particle.avgAttInSc(:,1:nLastDeployDepthImagSysLocal,2:config.nParticleAttributes,:,iRun) =...
            avgAttInScSeasonal(:,:,(2:config.nParticleAttributes),:); % attributes (the 1st attribute is discarded as it's the pNum transformed) 
        output.seasonal.particle.avgAttInTy(:,1:nLastDeployDepthImagSysLocal,1,:,iRun) =...
            pNumInTySeasonal; % particle number per Ty
        output.seasonal.particle.avgAttInTy(:,1:nLastDeployDepthImagSysLocal,2:5,:,iRun) =...
            avgAttInTySeasonal; % attributes
        if ~choiceSaveDiskSpace
            output.seasonal.particle.avgAttInVc(:,1:nLastDeployDepthImagSysLocal,1,:,iRun) =...
                pNumInVcSeasonal; % particle number per Vc
            output.seasonal.particle.avgAttInVc(:,1:nLastDeployDepthImagSysLocal,2:config.nParticleAttributes,:,iRun) =...
                avgAttInVcSeasonal(:,:,(2:config.nParticleAttributes),:); % attributes (the 1st attribute is discarded as it's the pNum transformed) 
        end  
        output.seasonal.fracMassAssocToTy(:,1:nLastDeployDepthImagSysLocal,:,iRun) = fracMassAssocToTySeasonal;
        output.seasonal.fracPocAssocToTy(:,1:nLastDeployDepthImagSysLocal,:,iRun) = fracPocAssocToTySeasonal; 

    end 

    % Process annual
    avgAttInScAnnual = aggregateParticleData('annual',avgAttInScMonthly,pNumInScMonthly,config);
    avgAttInTyAnnual = aggregateParticleData('annual',avgAttInTyMonthly,pNumInTyMonthly,config);
    if ~choiceSaveDiskSpace
        avgAttInVcAnnual = aggregateParticleData('annual',avgAttInVcMonthly,...
            pNumInVcMonthly,config);
    end
    fracMassAssocToTyAnnual = aggregateParticleData('annual',fracMassAssocToTyMonthly,...
        pNumInTyMonthly,config);
    fracPocAssocToTyAnnual = aggregateParticleData('annual',fracPocAssocToTyMonthly,...
        pNumInTyMonthly,config);

    % Handle NaN
    avgAttInScAnnual(isnan(avgAttInScAnnual)) = 0;
    avgAttInTyAnnual(isnan(avgAttInTyAnnual)) = 0;
    if ~choiceSaveDiskSpace
        avgAttInVcAnnual(isnan(avgAttInVcAnnual)) = 0;
    end
    fracMassAssocToTyAnnual(isnan(fracMassAssocToTyAnnual)) = 0;
    fracPocAssocToTyAnnual(isnan(fracPocAssocToTyAnnual)) = 0;

    % Save to output structure
    output.particle.avgAttInScAnnual(:,1:nLastDeployDepthImagSysLocal,1,iRun) =...
        pNumInScAnnual(:,:); % particle number per Sc 
    output.particle.avgAttInScAnnual(:,1:nLastDeployDepthImagSysLocal,(2:config.nParticleAttributes),iRun) =...
        avgAttInScAnnual(:,:,(2:config.nParticleAttributes)); % attributes (the 1st attribute is discarded as it's the pNum transformed) 
    output.particle.avgAttInTyAnnual(:,1:nLastDeployDepthImagSysLocal,1,iRun) =...
        pNumInTyAnnual(:,:); % particle number per Ty
    output.particle.avgAttInTyAnnual(:,1:nLastDeployDepthImagSysLocal,(2:5),iRun) =...
        avgAttInTyAnnual(:,:,:); % attributes
    if ~choiceSaveDiskSpace
        output.particle.avgAttInVcAnnual(:,1:nLastDeployDepthImagSysLocal,1,iRun) =...
            pNumInVcAnnual(:,:); % particle number per Vc
        output.particle.avgAttInVcAnnual(:,1:nLastDeployDepthImagSysLocal,(2:config.nParticleAttributes),iRun) =...
            avgAttInVcAnnual(:,:,(2:config.nParticleAttributes)); % attributes (the 1st attribute is discarded as it's the pNum transformed) 
    end
    output.particle.fracMassAssocToTyAnnual(:,1:nLastDeployDepthImagSysLocal,iRun) = fracMassAssocToTyAnnual(:,:); 
    output.particle.fracPocAssocToTyAnnual(:,1:nLastDeployDepthImagSysLocal,iRun) = fracPocAssocToTyAnnual(:,:);

    % .................................................................

    % Size and velocity stats

    if ~choiceSaveDiskSpace
        if choiceSeasonalOutput
            for iSeason = 1:config.nSeasons
                for iImagDepth = 1:nLastDeployDepthImagSysLocal
                    nonZeroSc = find(pNumInScSeasonal(:,iImagDepth,iSeason) ~= 0);
                    nonZeroVc = find(pNumInVcSeasonal(:,iImagDepth,iSeason) ~= 0);

                    pNumLargeParticlesAtDh = sum(pNumInScSeasonal(config.idxLargeParticles,iImagDepth,iSeason));
                    pNumFastParticlesAtDh = sum(pNumInVcSeasonal(config.idxFastParticles,iImagDepth,iSeason));
                    pNumAllParticlesAtDh = sum(pNumInScSeasonal(:,iImagDepth,iSeason)); % it's the same as sum(pNumInVcSeasonal(:,iId,iSeason))

                    output.seasonal.particle.sizeAndVeloStats(1,iImagDepth,iSeason,iRun) = pNumAllParticlesAtDh; % total number of particles at that depth horizon at that season
                    output.seasonal.particle.sizeAndVeloStats(2,iImagDepth,iSeason,iRun) = config.particleDiameterClasses(nonZeroSc(end)); % largest particle size class 
                    output.seasonal.particle.sizeAndVeloStats(3,iImagDepth,iSeason,iRun) = config.particleVelocityClasses(nonZeroVc(end)); % largest particle velocity class
                    output.seasonal.particle.sizeAndVeloStats(4,iImagDepth,iSeason,iRun) = pNumLargeParticlesAtDh/pNumAllParticlesAtDh; % fraction of large particles
                    output.seasonal.particle.sizeAndVeloStats(5,iImagDepth,iSeason,iRun) = pNumFastParticlesAtDh/pNumAllParticlesAtDh; % fraction of fast particles
                end
            end
        end

        for iImagDepth = 1:nLastDeployDepthImagSysLocal
            nonZeroSc = find(pNumInScAnnual(:,iImagDepth) ~= 0);            
            nonZeroVc = find(pNumInVcAnnual(:,iImagDepth) ~= 0); 

            pNumLargeParticlesAtDh = sum(pNumInScAnnual(config.idxLargeParticles,iImagDepth));            
            pNumFastParticlesAtDh = sum(pNumInVcAnnual(config.idxFastParticles,iImagDepth));
            pNumAllParticlesAtDh = sum(pNumInScAnnual(:,iImagDepth)); % it's the same as sum(pNumInVcAnnual(:,iId))

            output.particle.sizeAndVeloStatsAnnual(1,iImagDepth,iRun) = pNumAllParticlesAtDh; % total number of particles at that depth horizon
            output.particle.sizeAndVeloStatsAnnual(2,iImagDepth,iRun) = config.particleDiameterClasses(nonZeroSc(end)); % largest particle size class
            output.particle.sizeAndVeloStatsAnnual(3,iImagDepth,iRun) = config.particleVelocityClasses(nonZeroVc(end)); % largest particle velocity class
            output.particle.sizeAndVeloStatsAnnual(4,iImagDepth,iRun) = pNumLargeParticlesAtDh/pNumAllParticlesAtDh; % fraction of large particles
            output.particle.sizeAndVeloStatsAnnual(5,iImagDepth,iRun) = pNumFastParticlesAtDh/pNumAllParticlesAtDh; % fraction of fast particles
        end

    end
    
end % processEquilibriumFluxAndParticleAttributes

% *************************************************************************

function output = processSourcesMinusSinksAndAuxilliaryTerms(output,iRun,...
    config,modelSms,modelSmsIntegrated,modelAux,nWriteStepsEulerian,nLocalDepthLayers,...
    choiceSeasonalOutput,choiceSaveDiskSpace)
        
    % Initialise arrays to store data for the last 12 months of the 
    % simulation (equilibrium)
    smsMonthly = zeros(config.maxNoSMSterms,nLocalDepthLayers,12);
    auxMonthly = zeros(config.maxNoAuxTerms,nLocalDepthLayers,12);
    if ~choiceSaveDiskSpace
        smsMonthlyIntegrated = zeros(config.maxNoSMSterms,12);
    end
    for iMonth = 1:12           
        iSlice = (nWriteStepsEulerian-12) + iMonth;  
        smsMonthly(:,:,iMonth) = modelSms(:,:,iSlice); % mol m-3 d-1 (for a typical day of a particular month)
        auxMonthly(:,:,iMonth) = modelAux(:,:,iSlice);  
        if ~choiceSaveDiskSpace
            smsMonthlyIntegrated(:,iMonth) = modelSmsIntegrated(:,iSlice); % mol m-2 d-1 (for a typical day of a particular month)
        end
    end

    auxSeasonal = zeros(config.maxNoAuxTerms,nLocalDepthLayers,4);
    for iSeason = 1:config.nSeasons
        monthsInSeason = (3 * (iSeason - 1) + 1):(3 * iSeason); % get months for the season
        auxSeasonal(:,:,iSeason) = computeWeightedAverage(...
            auxMonthly(:,:,monthsInSeason),config.nDaysPerMonth(monthsInSeason),3);
    end

    % Seasonal calculation
    if choiceSeasonalOutput 
        smsSeasonal = zeros(config.maxNoSMSterms,nLocalDepthLayers,4);
        %auxSeasonal = zeros(config.maxNoAuxTerms,nLocalDepthLayers,4);
        for iSeason = 1:config.nSeasons
            monthsInSeason = (3 * (iSeason - 1) + 1):(3 * iSeason); % get months for the season
            smsSeasonal(:,:,iSeason) = computeWeightedAverage(...
                smsMonthly(:,:,monthsInSeason),config.nDaysPerMonth(monthsInSeason),3);
            % auxSeasonal(:,:,iSeason) = computeWeightedAverage(...
            %     auxMonthly(:,:,monthsInSeason),config.nDaysPerMonth(monthsInSeason),3);
        end
    end

    % Annual calculation
    smsAnnual = computeWeightedAverage(smsMonthly,config.nDaysPerMonth,3); 
    auxAnnual = computeWeightedAverage(auxMonthly,config.nDaysPerMonth,3);

    % Save results to output structure
    depthRange = 1:min(nLocalDepthLayers,config.maxNoDepthLayers); % determine the depth range to save
    output.sms.monthly(:,depthRange,:,iRun) = smsMonthly(:,depthRange, :);
    output.sms.annual(:,depthRange,iRun) = smsAnnual(:,depthRange);
    output.aux.annual(:,depthRange,iRun) = auxAnnual(:,depthRange);
    output.aux.monthly(:,depthRange,:,iRun) = auxMonthly(:,depthRange,:);
    output.aux.seasonal(:,depthRange,:,iRun) = auxSeasonal(:,depthRange,:);
    if choiceSeasonalOutput
        output.seasonal.sms(:,depthRange,:,iRun) = smsSeasonal(:,depthRange,:);
        % output.aux.seasonal(:,depthRange,:,iRun) = auxSeasonal(:,depthRange,:);
    end
    if ~choiceSaveDiskSpace
        output.sms.monthlyIntegrated(:,:,iRun) = smsMonthlyIntegrated;
    end

end % processSourcesMinusSinksAndAuxilliaryTerms

% ************************************************************************* 

function output = processSeafloorFluxAverages(output,iRun,config,...
    nLastDeployDepthSedTrapLocal,choiceSeasonalOutput,choiceSaveDiskSpace)
  
    fluxMonthly = output.flux.monthly(:,:,:,iRun);

    % Seasonal

    if choiceSeasonalOutput
        % I STILL HAVE TO CREATE IT!!!
    end

    % Annual

    if ~choiceSaveDiskSpace 

        % Initialise monthly fluxes and fractions
        nComponents = 3; % orgC, calcite, opal
        fluxSeabedMonthly = zeros(12, nComponents); % store fluxes for 12 months and 3 components
        fracSeabedMonthly = zeros(12, nComponents); % store fractions for 12 months and 3 components

        % Compute monthly fluxes and fractions
        for iMonth = 1:12
            fluxSeabedMonthly(iMonth,1) = sum(fluxMonthly(1:2,nLastDeployDepthSedTrapLocal,iMonth)); % POC + TEP
            fluxSeabedMonthly(iMonth,2) = fluxMonthly(3,nLastDeployDepthSedTrapLocal,iMonth); % Calcite
            fluxSeabedMonthly(iMonth,3) = fluxMonthly(4,nLastDeployDepthSedTrapLocal,iMonth); % Opal

            totFluxSeabedMonthly = sum(fluxSeabedMonthly(iMonth, :)); % total flux for the month
            fracSeabedMonthly(iMonth,:) = fluxSeabedMonthly(iMonth,:)./totFluxSeabedMonthly; % fractions
        end

        % Weighted annual sums
        weightedDays = config.nDaysPerMonth(:); % ensure column vector
        annualFlux = sum(fluxSeabedMonthly .* weightedDays, 1) / sum(weightedDays); 
        annualFrac = sum(fracSeabedMonthly .* weightedDays, 1) / sum(weightedDays); 

        % Compute standard deviations (weighted)
        stdFlux = std(fluxSeabedMonthly, weightedDays); 
        stdFrac = std(fracSeabedMonthly, weightedDays);

        % Save to files
        output.flux.seabedAnnual(1:3,1,1,iRun) = annualFlux; % annual fluxes
        output.flux.seabedAnnual(1:3,2,1,iRun) = annualFrac; % annual fractions
        output.flux.seabedAnnual(1:3,1,2,iRun) = stdFlux;    % flux standard deviations
        output.flux.seabedAnnual(1:3,2,2,iRun) = stdFrac;    % fraction standard deviations
    end 
        
end % processSeafloorFluxAverages
        
% *************************************************************************        

function output = calculateAverageDataForAllRunYears(output,iRun,config,...
    modelFlux,modelAvgAttSizeSpec,modelAvgAttSizeSpecSf,nLastDeployDepthSedTrapLocal,...
    nLastDeployDepthImagSysLocal)

    evolFluxMonthly = zeros(5,nLastDeployDepthSedTrapLocal,12,config.nYears);
    evolPnumInScMonthly = zeros(config.nDiameterClasses,nLastDeployDepthImagSysLocal,12,config.nYears);
    evolAvgAttInScMonthly = zeros(config.nDiameterClasses,nLastDeployDepthImagSysLocal,config.nParticleAttributes,12,config.nYears);

    % Group data from all years monthly
    iSlice = 0;
    for iYear = 1:config.nYears
        for iMonth = 1:12
            iSlice = iSlice + 1;

            % Flux data
            evolFluxMonthly(:,:,iMonth,iYear) = modelFlux(:,:,iSlice);

            % Particle number and attributes (water column and seafloor)
            evolPnumInScMonthly(:,1:nLocalImagSysDeployDepths,iMonth,iYear)      = modelAvgAttSizeSpec(:,:,1,iSlice);
            evolAvgAttInScMonthly(:,1:nLocalImagSysDeployDepths,:,iMonth,iYear)  = modelAvgAttSizeSpec(:,:,:,iSlice);
            evolPnumInScMonthly(:,nLastDeployDepthImagSysLocal,iMonth,iYear)     = modelAvgAttSizeSpecSf(:,1,1,iSlice);
            evolAvgAttInScMonthly(:,nLastDeployDepthImagSysLocal,:,iMonth,iYear) = modelAvgAttSizeSpecSf(:,1,:,iSlice);
        end
    end

    % Seasonal calculations
    evolFluxSeasonal = computeSeasonalData(evolFluxMonthly,config);
    evolPnumInScSeasonal = computeSeasonalData(evolPnumInScMonthly,config);
    evolAvgAttInScSeasonal = zeros(config.nDiameterClasses,nLastDeployDepthImagSysLocal,config.nParticleAttributes,4,config.nYears);
    for iSeason = 1:config.nSeasons
        monthsInSeason = (3 * (iSeason - 1) + 1):(3 * iSeason);
        for iYear = 1:config.nYears
            evolAvgAttInScSeasonal(:,:,:,iSeason,iYear) = sum(...
                evolAvgAttInScMonthly(:,:,:,monthsInSeason,iYear) .* ...
                evolPnumInScMonthly(:,:,monthsInSeason,iYear) .* ...
                reshape(nDaysPerMonth(monthsInSeason),1,1,[]), 3)...
                ./ evolPnumInScSeasonal(:,:,iSeason,iYear);
        end
    end

    % Handle NaN values
    evolPnumInScSeasonal(isnan(evolPnumInScSeasonal)) = 0;
    evolAvgAttInScSeasonal(isnan(evolAvgAttInScSeasonal)) = 0;

    % Annual calculations
    evolFluxAnnual = computeAnnualData(evolFluxMonthly,config);
    evolPnumInScAnnual = computeAnnualData(evolPnumInScMonthly,config);
    evolAvgAttInScAnnual = zeros(config.nDiameterClasses,nLastDeployDepthImagSysLocal,config.nParticleAttributes,config.nYears);
    for iYear = 1:config.nYears
        evolAvgAttInScAnnual(:,:,:,iYear) = sum(...
            evolAvgAttInScMonthly(:,:,:,:,iYear) .* ...
            evolPnumInScMonthly(:,:,:,iYear) .* ...
            reshape(nDaysPerMonth,1,1,[]), 3)...
            ./ evolPnumInScAnnual(:,:,iYear);
    end

    % Handle NaN values
    evolFluxAnnual(isnan(evolFluxAnnual)) = 0;
    evolPnumInScAnnual(isnan(evolPnumInScAnnual)) = 0;
    evolAvgAttInScAnnual(isnan(evolAvgAttInScAnnual)) = 0;

    % Save to output structure
    output.evolution.flux.seasonal(:,1:nLastDeployDepthSedTrapLocal,:,:,iRun)                 = evolFluxSeasonal;
    output.evolution.flux.annual(:,1:nLastDeployDepthSedTrapLocal,:,iRun)                     = evolFluxAnnual; 
    output.evolution.particle.avgAttInScSeasonal(:,1:nLastDeployDepthImagSysLocal,:,:,:,iRun) = evolAvgAttInScSeasonal(:,:,(2:config.nParticleAttributes),:,:);
    output.evolution.particle.avgAttInScAnnual(:,1:nLastDeployDepthImagSysLocal,:,:,iRun)     = evolAvgAttInScAnnual(:,:,(2:config.nParticleAttributes),:);
    output.evolution.particle.numberInScAnnual(:,1:nLastDeployDepthImagSysLocal,:,iRun)       = evolPnumInScAnnual;
            
end % calculateAverageDataForAllRunYears
            
% *************************************************************************

function output = processParticleDataBySizeAndVelocityClass(output,iRun,...
    config,nLastDeployDepthImagSysLocal,choiceSeasonalOutput,choiceSaveDiskSpace)

    pNumInScMonthly = output.particle.avgAttInScMonthly(:,:,1,:,iRun);
    pNumInVcMonthly = output.particle.avgAttInVcMonthly(:,:,1,:,iRun);
    avgAttInScMonthly = output.particle.avgAttInScMonthly(:,:,2:config.nParticleAttributes,:,iRun);
    avgAttInVcMonthly = output.particle.avgAttInVcMonthly(:,:,2:config.nParticleAttributes,:,iRun);
       
    % if choiceSeasonalOutput

        % 1st dim (class): 1=all Sc, 2=small Sc, 3=large Sc, 4=all Vc, 5=slow Vc, 6=fast Vc
        totPnumAtDhMonthlyAdd    = zeros(6,nLastDeployDepthImagSysLocal,4); 
        totPnumAtDhSeasonal      = zeros(6,nLastDeployDepthImagSysLocal,4);
        % 1st dim (class): 1=all Sc, 2=small Sc, 3=large Sc
        avgAttInScAtDhMonthlyAdd = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1,4); % 4th dim: (1:config.nParticleAttributes) for atts, (+1) to accommodate average velocity/size
        avgAttInScAtDhSeasonal   = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1,4);
        % 1st dim (class): 1=all Vc, 2=slow Vc, 3=fast Vc
        avgAttInVcAtDhMonthlyAdd = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1,4);
        avgAttInVcAtDhSeasonal   = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1,4);

        % Seasonal calculation loop
        iSeason = 1;
        for iMonth = 1:12
            for iImagDepth = 1:nLastDeployDepthImagSysLocal

                % Sum particle counts across classes
                for iClass = 1:6
                    classRange = getClassRange(iClass,config);
                    if iClass <= 3
                        % Small/large particles (Sc)
                        pNumInClassRange = sum(pNumInScMonthly(classRange,iImagDepth,iMonth));
                    else
                        % Slow/fast velocities (Vc)
                        pNumInClassRange = sum(pNumInVcMonthly(classRange,iImagDepth,iMonth));
                    end 
                    totPnumAtDhMonthlyAdd(iClass,iImagDepth,iSeason) =...
                        totPnumAtDhMonthlyAdd(iClass,iImagDepth,iSeason) + (pNumInClassRange*config.nDaysPerMonth(iMonth));
                end

                if sum(pNumInScMonthly(:,iImagDepth,iMonth)) > 0 % if there are particles                 

                    % Average attributes calculation
                    for iAtt = 1:(config.nParticleAttributes-1) 
                        for iClass = 1:6
                            classRange = getClassRange(iClass,config);
                            if iClass <= 3
                                % Small/large particles (Sc)
                                avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,iAtt,iSeason) =... 
                                    avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,iAtt,iSeason)...
                                    + sum( avgAttInScMonthly(classRange,iImagDepth,iAtt,iMonth)...
                                        .* pNumInScMonthly(classRange,iImagDepth,iMonth).*config.nDaysPerMonth(iMonth) );
                            else
                                % Slow/fast velocities (Vc)
                                avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,iAtt,iSeason) =...
                                    avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,iAtt,iSeason)...
                                    + sum( avgAttInVcMonthly(classRange,iImagDepth,iAtt,iMonth)...
                                        .* pNumInVcMonthly(classRange,iImagDepth,iMonth).*config.nDaysPerMonth(iMonth) );
                            end 

                        end
                    end

                    % Special case for size and velocity attributes (compute average diameter/velocity)
                    for iClass = 1:6
                        classRange = getClassRange(iClass,config);
                        if iClass <= 3
                            avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,config.nParticleAttributes,iSeason) =...
                                avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,config.nParticleAttributes,iSeason)...
                                + sum( config.particleDiameterClasses(classRange)...
                                    .* pNumInScMonthly(classRange,iImagDepth,iMonth)*config.nDaysPerMonth(iMonth) );
                        else
                            avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,config.nParticleAttributes,iSeason) =...
                                avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,config.nParticleAttributes,iSeason)...
                                + sum( config.particleVelocityClasses(classRange)...
                                    .* pNumInVcMonthly(classRange,iImagDepth,iMonth).*config.nDaysPerMonth(iMonth) ); 
                        end
                    end

                end % checking if there are particles
            end % iImagDepth

            % Update season every 3 months
            if mod(iMonth,3) == 0
                iSeason = iSeason + 1;
            end 

        end % iMonth

        % Average over seasonal data
        for iSeason = 1:config.nSeasons

            % In order to average the total number of particles in each big
            % class, we divide the monthly pondered adddition by the sum of
            % the number of days they have occurred, e.g. (1000*31 + 800*28
            % + 950*30)/(31+28+30)
            totPnumAtDhSeasonal(:,:,iSeason) = totPnumAtDhMonthlyAdd(:,:,iSeason)./config.nDaysPerSeason(iSeason);

            % In order to verage the attributes, we divide the pondered
            % addition by the number of particles (i.e., 2.8*1000 + 1.7*100
            % + 0.8*10000)/(1000+100+10000))
            for iAtt = 1:config.nParticleAttributes
                avgAttInScAtDhSeasonal(:,:,iAtt,iSeason) = avgAttInScAtDhMonthlyAdd(:,:,iAtt,iSeason)./totPnumAtDhMonthlyAdd((1:3),:,iSeason);
                avgAttInVcAtDhSeasonal(:,:,iAtt,iSeason) = avgAttInVcAtDhMonthlyAdd(:,:,iAtt,iSeason)./totPnumAtDhMonthlyAdd((4:6),:,iSeason);
            end

        end % iSeason

        % Handle NaN values
        totPnumAtDhSeasonal(isnan(totPnumAtDhSeasonal)) = 0;
        avgAttInScAtDhSeasonal(isnan(avgAttInScAtDhSeasonal)) = 0;
        avgAttInVcAtDhSeasonal(isnan(avgAttInVcAtDhSeasonal)) = 0;

        % Save to output structure
        output.seasonal.particle.avgAttScAtDepth(:,1:nLastDeployDepthImagSysLocal,1,:,iRun) ...
            = totPnumAtDhSeasonal((1:3),:,:); % total number of particles across size classes by size categories (all, small, large)
        output.seasonal.particle.avgAttVcAtDepth(:,1:nLastDeployDepthImagSysLocal,1,:,iRun) ...
             = totPnumAtDhSeasonal((4:6),:,:); % total number of particles across velocity classes by velocity categories (all, slow, fast)
        
        output.seasonal.particle.avgAttScAtDepth(:,1:nLastDeployDepthImagSysLocal,(2:config.nParticleAttributes+1),:,iRun) ...
            = avgAttInScAtDhSeasonal(:,:,(1:config.nParticleAttributes),:); % attributes across size classes by size categories (all, small, large)
        output.seasonal.particle.avgAttVcAtDepth(:,1:nLastDeployDepthImagSysLocal,(2:config.nParticleAttributes+1),:,iRun) ...
             = avgAttInVcAtDhSeasonal(:,:,(1:config.nParticleAttributes),:); % attributes velocity classes by velocity categories (all, slow, fast)
   

    % end % choiceSeasonalOutput
    
    % .....................................................................

    % Calculate annual values

    % 1st dim (class): 1=all Sc, 2=small Sc, 3=large Sc, 4=all Vc, 5=slow Vc, 6=fast Vc
    totPnumAtDhMonthlyAdd    = zeros(6,nLastDeployDepthImagSysLocal); 
    totalPnumAtDhAnnual      = zeros(6,nLastDeployDepthImagSysLocal);
    % 1st dim (class): 1=all Sc, 2=small Sc, 3=large Sc
    avgAttInScAtDhMonthlyAdd = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1); % 3rd dim: (1:config.nParticleAttributes) for atts, (+1) to accommodate average velocity/size
    avgAttInScAtDhAnnual     = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1);
    % 1st dim (class): 1=all Vc, 2=slow Vc, 3=fast Vc
    avgAttInVcAtDhMonthlyAdd = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1);
    avgAttInVcAtDhAnnual     = zeros(3,nLastDeployDepthImagSysLocal,config.nParticleAttributes+1);

    % Add the values of all months together
    for iMonth = 1:12 
        for iImagDepth = 1:nLastDeployDepthImagSysLocal

            % Sum particle counts across classes
            for iClass = 1:6
                classRange = getClassRange(iClass,config);
                if iClass <= 3
                    % Small/large particles (Sc)
                    pNumInClassRange = sum(pNumInScMonthly(classRange,iImagDepth,iMonth));
                else
                    % Slow/fast velocities (Vc)
                    pNumInClassRange = sum(pNumInVcMonthly(classRange,iImagDepth,iMonth));
                end   
                totPnumAtDhMonthlyAdd(iClass,iImagDepth) =...
                    totPnumAtDhMonthlyAdd(iClass,iImagDepth) + (pNumInClassRange.*config.nDaysPerMonth(iMonth));
            end

            if sum(pNumInScMonthly(:,iImagDepth,iMonth)) > 0 % if there are particles               
                for iAtt = 1:(config.nParticleAttributes-1) 
                    for iClass = 1:6
                        classRange = getClassRange(iClass,config);
                        if iClass <= 3
                            % Small/large particles (Sc)              
                            avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,iAtt) =...
                                avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,iAtt)...
                                    + sum( avgAttInScMonthly(classRange,iImagDepth,iAtt,iMonth)...
                                        .*pNumInScMonthly(classRange,iImagDepth,iMonth).*config.nDaysPerMonth(iMonth) ); 
                        else
                            % Slow/fast velocities (Vc)
                            avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,iAtt) =... 
                                avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,iAtt)...
                                    + sum( avgAttInVcMonthly(classRange,iImagDepth,iAtt,iMonth)...
                                        .*pNumInVcMonthly(classRange,iImagDepth,iMonth).*config.nDaysPerMonth(iMonth) );
                        end
                    end
                end % iAtt

                % Special case for size and velocity attributes (compute average diameter/velocity)
                for iClass = 1:6
                    classRange = getClassRange(iClass,config);
                    if iClass <= 3
                        avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,config.nParticleAttributes) ...
                            = avgAttInScAtDhMonthlyAdd(iClass,iImagDepth,config.nParticleAttributes)...
                                + sum( config.particleDiameterClasses(classRange)...
                                    .*pNumInScMonthly(classRange,iImagDepth,iMonth).*config.nDaysPerMonth(iMonth) );
                    else  
                        avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,config.nParticleAttributes) ...
                            = avgAttInVcAtDhMonthlyAdd(iClass-3,iImagDepth,config.nParticleAttributes)...
                                + sum( config.particleVelocityClasses(classRange)...
                                    .*pNumInVcMonthly(classRange,iImagDepth,iMonth).*config.nDaysPerMonth(iMonth) );
                    end
                end
            end % if there are particles
        end % iImagDepth    
    end % iMonth

    % In order to average the total number of particles in each big
    % class, we divide the monthly pondered adddition by the sum of
    % the number of days they have occurred, e.g. (1000*31 + 800*28
    % + 950*30)/(31+28+30)
    totalPnumAtDhAnnual(:,:) = totPnumAtDhMonthlyAdd(:,:)./sum(config.nDaysPerMonth(:));

    % Average attribute per big size class and velocity class
    for iAtt = 1:config.nParticleAttributes
       avgAttInScAtDhAnnual(:,:,iAtt) = avgAttInScAtDhMonthlyAdd(:,:,iAtt)./totPnumAtDhMonthlyAdd((1:3),:);
       avgAttInVcAtDhAnnual(:,:,iAtt) = avgAttInVcAtDhMonthlyAdd(:,:,iAtt)./totPnumAtDhMonthlyAdd((4:6),:);
    end

    % Handle NaN values
    totalPnumAtDhAnnual(isnan(totalPnumAtDhAnnual)) = 0;
    avgAttInScAtDhAnnual(isnan(avgAttInScAtDhAnnual)) = 0;
    avgAttInVcAtDhAnnual(isnan(avgAttInVcAtDhAnnual)) = 0;

    % Save to output structure
    output.particle.avgAttScAtDepthAnnual(:,1:nLastDeployDepthImagSysLocal,1,iRun) ...
        = totalPnumAtDhAnnual((1:3),:); % total number of particles across size classes by size categories (all, small, large)
    output.particle.avgAttScAtDepthAnnual(:,1:nLastDeployDepthImagSysLocal,(2:config.nParticleAttributes+1),iRun) ...
        = avgAttInScAtDhAnnual(:,:,(1:config.nParticleAttributes)); % attributes across size classes by size categories (all, small, large)
    
    if ~choiceSaveDiskSpace
        output.particle.avgAttVcAtDepthAnnual(:,1:nLastDeployDepthImagSysLocal,1,iRun) ...
            = totalPnumAtDhAnnual((4:6),:); % total number of particles across velocity classes by velocity categories (all, slow, fast)
        output.particle.avgAttVcAtDepthAnnual(:,1:nLastDeployDepthImagSysLocal,(2:config.nParticleAttributes+1),iRun) ...
            = avgAttInVcAtDhAnnual(:,:,(2:config.nParticleAttributes)); % attributes across velocity classes by velocity categories (all, slow, fast)
    end

end % processParticleDataBySizeAndVelocityClass

% *************************************************************************

end % readSlamsOutput