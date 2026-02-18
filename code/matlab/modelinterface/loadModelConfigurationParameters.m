function config = loadModelConfigurationParameters(fullpathModelInputDataDir,...
    filenameSlamsRunGrid,filenameSlamsNumDepthLayers,choiceTypeGridDomain)

% LOADMODELCONFIGURATIONPARAMETERS Structure that contains fixed quantities
% used to pre-process, run and post-process SLAMS. Make sure some sections
% in this script are up-to-date with the following configuration files:
%   modelparameters.F90
%   namelist
%   modeloutput.F90

config = struct();

%% Time and depth stepping (in namelist, &ConfigurationParameters)

config.nYears = 6;
config.nTimestepsPerDay = 3; % 28800 s = 8 h = 1/3 of a day
config.maxSimulatedDepth = 2000; % from nDepthLayers = 400
config.depthLayerArea = 0.10; % m2
config.defaultDepthLayerThickness = 10; % m
config.depthLayerVolume = config.depthLayerArea * config.defaultDepthLayerThickness;
config.nTimestepsPerYear = config.nTimestepsPerDay * 365;

%% Deployment depths for collecting fluxes and imaging particles (in modelparameters.F90, GetDepthsForParticleDataCollection)

[config.availSedTrapDeployDepths,config.maxNoSedTrapDeployDepths] = defineDeploymentDepths();
[config.availImagSysDeployDepths,config.maxNoImagSysDeployDepths] = defineDeploymentDepths();

%% Grid-related parameters (read from external files)

gridFile = fullfile(fullpathModelInputDataDir,filenameSlamsRunGrid);
nDepthLayersFile = fullfile(fullpathModelInputDataDir,filenameSlamsNumDepthLayers);

% Load number of depth layers
if isfile(nDepthLayersFile)
    nDepthLayersUsedID = fopen(nDepthLayersFile,'r');
    config.nDepthLayersUsed = fscanf(nDepthLayersUsedID,'%d\n'); 
    fclose(nDepthLayersUsedID);
else
    fprintf('Depth file "%s" does not exist yet... creating it.\n', nDepthLayersFile);
end

% Load grid
if isfile(gridFile)
    if choiceTypeGridDomain == 1 || choiceTypeGridDomain == 3 % global
        load(gridFile,'Xbb','Ybb','Zsb','x','y','ixBb','iyBb')
    elseif choiceTypeGridDomain == 2 % local
        load(gridFile,'Xbb','Ybb','Zsb')
    end

    config.nLocs = length(Ybb);
    config.gridLats = Ybb;
    config.gridLons = Xbb;
    gridDepths = Zsb;
    
    if choiceTypeGridDomain == 1 || choiceTypeGridDomain == 3 % global
        config.nLats = length(y);
        config.nLons = length(x);
        config.lats = y;
        config.lons = x;
        config.ixBb = ixBb;
        config.iyBb = iyBb;
        config.lonIncrement = x(2) - x(1); % degrees longitude worth one cell
        config.latIncrement = y(2) - y(1); % degrees latitude worth one cell

        % Verify uniformity
        lonDiffs = diff(x); % differences between consecutive longitude values
        latDiffs = diff(y); % differences between consecutive latitude values
        config.isLonUniform = all(abs(lonDiffs - config.lonIncrement) < 1e-10); % allow small tolerance
        config.isLatUniform = all(abs(latDiffs - config.lonIncrement) < 1e-10);

    else % local
        config.lonIncrement = 1; % degrees longitude 
        config.latIncrement = 1; % degrees latitude
        config.nLats = length(Ybb);
        config.nLons = length(Xbb);
        config.lats = Ybb;
        config.lons = Xbb;
    end

    % Largest number of steps that fit entirely within the maximum boundary 
    % without exceeding it
    maxBoundary = 20.0; % maximum boundary in degrees
    config.maxIterForNeighbourhoodExpansion = floor(maxBoundary/config.lonIncrement); % number of full steps
     
else
    error('Grid file "%s" does not exist.', filenameSlamsRunGrid);
    return; % exit the script or function 
end

% Mask depths greater than maxSimulatedDepth for performance improvement
gridDepths(gridDepths > config.maxSimulatedDepth) = NaN;
config.gridDepths = gridDepths;

% Compute maximum number of valid depth layers
config.maxNoDepthLayers = max(sum(~isnan(gridDepths), 1));

% Rearrange the depth grid from sequence of values to latitude x longitude
% structure
if choiceTypeGridDomain == 1 || choiceTypeGridDomain == 3 % global
    config.geoDepths = NaN(config.nLats,config.nLons,config.maxNoDepthLayers);
    for iLoc = 1:config.nLocs
        iLat = config.iyBb(iLoc);
        iLon = config.ixBb(iLoc);
        validDepths = config.gridDepths(~isnan(config.gridDepths(:,iLoc)),iLoc);
        config.geoDepths(iLat,iLon,1:numel(validDepths)) = validDepths; 
    end
end

%% Particle attributes (in modelparameters.F90)

% Categories defined in pstruct.h
config.nParticleTypes = 9;
% (1) faecal                   (2) non-faecal aggregate
% (3) living phytoplankton     (4) dead phytoplankton
% (5) TEP                      (6) clay 
% (7) zooplankton carcass      (8) inorganic, single, non-faecal
% (9) organic, single, non-faecal that is not dead phyto (fragment)

% Categories defined in GetParticleAverageAttributes in modeloutput.F90
config.nParticleAttributes = 14;
% (1) no. particles            (2) density
% (3) size/velocity            (4) stickiness 
% (5) porosity                 (6) solid volume 
% (7) fractal dimension        (8) radius PP
% (9) depth                    (10) moles orgC 
% (11) moles TEPC              (12) moles opal
% (13) moles calcite           (14) moles clay

% Categories defined in EstablishParticleSizeAndVeloCategories in modelparameters.F90
config.nDiameterClasses = 16;
config.nVelocityClasses = 16;

config.particleDiameterClasses = 1 .* 2.^(0:config.nDiameterClasses-1)'; % doubling sequence starting from 1 um
config.particleVelocityClasses = 0.1 .* 2.^(0:config.nVelocityClasses-1)'; % doubling sequence starting from 0.1 m d-1

% For post-processing output
config.diameterLargeParticle = 150; % 150 um, diameter (as in Giering et al., 2020)
config.velocityFastParticle  = 100; % 100 m d-1 (McDonnell and Buesseler (2010) suggest 150 m d-1)

config.idxLargeScThreshold = find(config.particleDiameterClasses <= config.diameterLargeParticle, 1, 'last');
config.idxFastVcThreshold = find(config.particleVelocityClasses <= config.velocityFastParticle, 1, 'last');

config.idxLargeParticles = find(config.particleDiameterClasses >= config.diameterLargeParticle);  
config.idxFastParticles = find(config.particleVelocityClasses >= config.velocityFastParticle);

% Other diameter-and-velocity-related attributes
config.particleDiameterClassBounds = [0.2; config.particleDiameterClasses];
config.particleDiameterClassMiddle = (config.particleDiameterClassBounds(2:end) + config.particleDiameterClassBounds(1:end-1)) ./ 2;
config.particleDiameterClassWidth = diff(config.particleDiameterClassBounds);

config.particleVelocityClassBounds = [0; config.particleVelocityClasses];
config.particleVelocityClassMiddle = (config.particleVelocityClassBounds(2:end) + config.particleVelocityClassBounds(1:end-1)) ./ 2;
config.particleVelocityClassWidth = diff(config.particleVelocityClassBounds);

%% Model PFT constants (in the namelist, &FixedParameters)

config.pft.q10 = 1.65 * ones(4,1); % all PFTs have the same Q10 value
config.pft.tempRef = 0;
config.pft.muMax0degPft = [0.70, 0.20, 0.50,  0.25]; % in this order: diatoms, dino, coccos, pico
config.pft.kNO3Pft      = [2.0,  3.0,  0.80,  0.25];
config.pft.kPO4Pft      = [0.20, 0.30, 0.005, 0.004];
config.pft.kSiDiat = 4.0;  % supported by Sarthou et al. 2005
config.pft.chlToCarbonMax = 2.5e-2;
config.pft.alphaChlSpecificPft = [2.0e-5, 0.7e-5, 0.3e-5, 0.7e-5];
config.pft.muThreshold = 2.0; % d-1

%% Indexes to SMS terms (mol C m-3 d-1) (in modelparameters.F90)

config.maxNoSMSterms = 50;

config.sms.idxPrimProdOrgC     = 1; % C taken up by phytoplankton
config.sms.idxPrimProdCaCO3    = 2; % CO3= taken up by coccos
config.sms.idxPrimProdOpal     = 3; % SiOH4 taken up by diatoms					  
config.sms.idxDepoClay	       = 4; % aeolian dust deposition
config.sms.idxProdTepPhyto     = 5; % TEC released by phytoplankton during photosynthesis
config.sms.idxProdTepMicrob    = 6; % TEC released by microbes
config.sms.idxPhotoTepC        = 7; % DOC released due to photolysis	  
config.sms.idxZooSolubOrgC     = 8; % DOC released due to zooplankton messy feeding
config.sms.idxZooSolubTepC     = 9; % dissolved TEC released due to zooplankton messy feeding
config.sms.idxZooSolubCaCO3    = 10; % CO3= released due to zooplankton messy feeding
config.sms.idxZooSolubOpal     = 11; % SiOH4 released due to zooplankton messy feeding
config.sms.idxZooSolubClay     = 12; % dissolved clay released due to zooplankton messy feeding
config.sms.idxZooIngestOrgC    = 13; % zooplankton ingestion of POC
config.sms.idxZooIngestTepC    = 14; % zooplankton ingestion of TEC
config.sms.idxZooRespOrgC      = 15; % C (CO2) released due to zooplankton respiration of POC
config.sms.idxZooRespTepC      = 16; % C (CO2) released due to zooplankton respiration of TEC
config.sms.idxZooEgestOrgC     = 17; % zooplankton egestion of POC
config.sms.idxZooEgestTepC     = 18; % zooplankton egestion of TEC
config.sms.idxZooDissolCaCO3   = 19; % CO3= release due to biotic dissolution of CaCO3
config.sms.idxZooExcretOrgC    = 20; % DOC release due to zooplankton excretion
config.sms.idxZooExcretTepC    = 21; % dissolved TEC release due to zooplankton excretion
config.sms.idxZooDeathOrgC     = 22; % POC released as zooplankton dead bodies
config.sms.idxZooDeathCaCO3    = 23; % CaCO3 released as zooplankton dead bodies
config.sms.idxZooDeathOpal     = 24; % opal released as zooplankton dead bodies
config.sms.idxMicrobRespOrgC   = 25; % C (CO2) released due to microbial respiration of POC
config.sms.idxMicrobRespTepC   = 26; % C (CO2) released due to microbial respiration of TEC
config.sms.idxMicrobSolubOrgC  = 27; % DOC released due to microbial solubilisation
config.sms.idxMicrobSolubTepC  = 28; % dissolved TEC released due to microbial solubilisation
config.sms.idxMicrobSolubCaCO3 = 29; % CO3= released due to abiotic solubilisation
config.sms.idxMicrobSolubOpal  = 30; % SiOH4 released due to abiotic solubilisation
config.sms.idxMicrobSolubClay  = 31; % dissolved clay released due to abiotic solubilisation					  
config.sms.idxDissolCaCO3      = 32; % CO3= released due to abiotic dissolution of CaCO3	
config.sms.idxDissolOpal       = 33; % SiOH4 released due to abiotic dissolution of opal

%% Indexes to auxilliary terms (in modelparameters.F90)

config.maxNoAuxTerms = 40;

config.aux.idxZooBiomass             	= 1; % mol C m-3
config.aux.idxZooNumber              	= 2; % ind. m-3
config.aux.idxZooDeadBiomass         	= 3; % mol C m-3
config.aux.idxDiatBiomass            	= 4; % mol C m-3
config.aux.idxFlagelBiomass             = 5; % mol C m-3
config.aux.idxCoccoBiomass              = 6; % mol C m-3
config.aux.idxPicoBiomass            	= 7; % mol C m-3
config.aux.idxFreshDiatCellQuota     	= 8; % mol C cell-1
config.aux.idxFreshFlagelCellQuota      = 9; % mol C cell-1
config.aux.idxFreshCoccoCellQuota    	= 10; % mol C cell-1
config.aux.idxFreshPicoCellQuota     	= 11; % mol C cell-1				  
config.aux.idxAvgProbDiat	            = 12;
config.aux.idxAvgProbFlagel	            = 13;
config.aux.idxAvgProbCocco	            = 14;
config.aux.idxAvgProbPico	            = 15;
config.aux.idxProfProbDiat              = 16;
config.aux.idxProfProbFlagel            = 17;
config.aux.idxProfProbCocco             = 18;
config.aux.idxProfProbPico              = 19;
config.aux.idxZooSpecRespRate        	= 20; % s-1
config.aux.idxMicrobSpecRespRate     	= 21; % s-1
config.aux.idxNightDvmUpperBound     	= 22; % m
config.aux.idxNightDvmLowerBound     	= 23; % m
config.aux.idxDayDvmUpperBound       	= 24; % m
config.aux.idxDayDvmLowerBound	        = 25; % m
config.aux.idxEuphoticDepth          	= 26; % m
config.aux.idxCollisionKernel        	= 27; % m3 s-1
config.aux.idxBrownianKernel	        = 28; % m3 s-1
config.aux.idxShearKernel	          	= 29; % m3 s-1
config.aux.idxSettlingKernel	      	= 30; % m3 s-1
config.aux.idxNumClusterPairsEvalCoagu  = 31; % num.
config.aux.idxNumPxCmore		        = 32;
config.aux.idxNumPxCless		        = 33;
config.aux.idxCoagulationSuccess     	= 34;
config.aux.idxCoagulationProbability 	= 35;
config.aux.idxNumParticlesEvalEncounter = 36;
config.aux.idxNumParticlesEncountered   = 37;	
config.aux.idxNumParticlesFragmentedZoo = 38;
config.aux.idxNumParticlesIngestedZoo   = 39;
config.aux.idxNumParticlesOmittedZoo    = 40;

%% Model input filenames (in modelparameters.F90)

% Model input forcing 
config.filenameModelPar0         = 'in_par0.bin';
config.filenameModelNpp          = 'in_npp.bin';
config.filenameModelChla         = 'in_chla.bin';
config.filenameModelRho          = 'in_rho.bin';
config.filenameModelOmegaCalc    = 'in_omegacalcite.bin';
config.filenameModelNit          = 'in_nitrate.bin'; 
config.filenameModelSil          = 'in_silicicacid.bin';
config.filenameModelPhos         = 'in_phosphate.bin';
config.filenameModelOxy          = 'in_oxygen.bin';
config.filenameModelTemp         = 'in_temperature.bin';
config.filenameModelDust         = 'in_dust.bin';
config.filenameModelMld          = 'in_mld.bin';
config.filenameModelMesozoo      = 'in_mesozooplankton.bin';
config.filenameModelDynVisco     = 'in_dynamicviscosity.bin';

% Upper boundary layers and lower boundary layers
config.filenameModelZub = 'in_depthupperbounds.bin';
config.filenameModelZlb = 'in_depthlowerbounds.bin';

%% Model output filenames (in modelparameters.F90)

% Control and particle stats from the final snapshot
config.filenameOutputStatus         = 'program_status.tmp';
config.filenameOutputControl        = 'out_control.bin';
config.filenameOutputNumClusters    = 'out_stats_nClusters.bin';
config.filenameOutputNumParticles   = 'out_stats_nParticles.bin';
config.filenameOutputClusterInteger = 'out_cluster_i.bin';
config.filenameOutputClusterReal    = 'out_cluster_r.bin';
config.filenameInstAvgAttSizeClass   = 'out_inst_avgparticleatt_sizeclass.bin';
config.filenameInstAvgAttVeloClass   = 'out_inst_avgparticleatt_veloclass.bin';
config.filenameInstAvgAttSizeClassSf = 'out_inst_avgparticleatt_sizeclass_seafloor.bin';
config.filenameInstAvgAttVeloClassSf = 'out_inst_avgparticleatt_veloclass_seafloor.bin';

% Fluxes
config.filenameOutputFlux   = 'out_flux.bin';
config.filenameOutputFluxSf = 'out_flux_seafloor.bin';

% Particle properties
config.filenameOutputAvgAttsSizeClass      = 'out_avgparticleatt_sizeclass.bin';
config.filenameOutputAvgAttsSizeClassSf    = 'out_avgparticleatt_sizeclass_seafloor.bin';
config.filenameOutputAvgAttsVeloClass      = 'out_avgparticleatt_veloclass.bin';
config.filenameOutputAvgAttsVeloClassSf    = 'out_avgparticleatt_veloclass_seafloor.bin';
config.filenameOutputAvgAttsParticleType   = 'out_avgparticleatt_maintype.bin';
config.filenameOutputAvgAttsParticleTypeSf = 'out_avgparticleatt_maintype_seafloor.bin';

% Sources-minus-sinks terms and auxilliary variables
config.filenameOutputSms           = 'out_sms.bin';
config.filenameOutputSmsIntegrated = 'out_sms_integrated.bin';
config.filenameOutputAux           = 'out_aux.bin';

%% Indexes to attribute array (in modeloutput.F90, WriteSnapshotsOfParticleAttributes)

config.cluster.id               = 1;
config.cluster.initType			= 2;
config.cluster.initPft		    = 3;
config.cluster.phase			= 4;
config.cluster.living			= 5;
config.cluster.faecal			= 6;
config.cluster.tstepCreat		= 7;
config.cluster.nPxC				= 8;
config.cluster.nPpxP			= 9;
config.cluster.nPpxC			= 10;
config.cluster.molesOrgC		= 11;
config.cluster.molesTepC		= 12;
config.cluster.molesOpal		= 13;
config.cluster.molesCalcite	    = 14;
config.cluster.molesClay		= 15;
config.cluster.massOrgMatter	= 16;
config.cluster.massTep			= 17;
config.cluster.depth			= 18;
config.cluster.mass				= 19;
config.cluster.solidVolume		= 20;
config.cluster.density			= 21;
config.cluster.radius			= 22;
config.cluster.radiusPp			= 23;
config.cluster.porosity			= 24;
config.cluster.stickiness		= 25;
config.cluster.velocity			= 26;

%% Other

config.maxNoClusters = 5e5; % in modelparameters.F90
config.fluxDetectionLimit = 1e-3; % mg m-2 d-1
config.maxNoClustersCutoff = 1e4;

%% Mapping properties

config.mappingProps = struct(...
    'mapLats', linspace(-90,90,config.nLats),... 
    'mapLons', linspace(-180,180,config.nLons),...
    'myColourMapOceanVars', brewermap(1000,'*RdYlBu'),... %[ones(1,3); jet(1000)]
    'myColourMapParticles', brewermap(1000,'*Spectral'),... % 'Spectral' is the equivalent of 'jet'
    'labelMonths', {{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'}},...
    'labelThreePfts', {{'Microphytoplankton','Nanophytoplankton','Picophytoplankton'}},...
    'labelFourPfts', {{'Diatoms','Large non-mineralising phytoplankton','Coccolithophores','Picophytoplankton'}}...
);

%% Molar mass and conversion factors

config.molarMassCarbon = 12.011; % g C / mol C
config.molarMassCaCO3 = 100.1;
config.molarMassBiogenicSilica = 67.3;
config.molarMassSilicon = 28.0855; % g mol-1
config.molarVolumeOxygen = 22.392; % L / mol (at STP)
config.dustToClayRatio = 0.5; % 0.5 g clay / 1 g dust (after Table 1 in Journet et al. 2008)

%% Depth horizon parameters (for post-processing calculations)

config.zmidmesoDef = 500;
config.zmesoDef = 1000;
config.zbathyDef = 2000;
config.errorFractionZeu = 0.10; % 10%, based on McKinna et al. 2019
config.errorZmeso = 20; % m
config.errorZbathy = 10; % m
config.randErrorFractionPocFlux = 0.30; % 30% (Buesseler et al. 2000, Buesseler et al. 2007, Stanley et al. 2004)
config.sysErrorFractionPocFlux  = 0.10; % 10%, based on a literature review
config.randErrorFractionPicFlux = 0.30; % 30% (Buesseler et al. 2000, Buesseler et al. 2007, Stanley et al. 2004)
config.sysErrorFractionPicFlux  = 0.10; % 10%, based on a literature review
config.randErrorFractionBSiFlux = 0.30; % 30% (Buesseler et al. 2000, Buesseler et al. 2007, Stanley et al. 2004)
config.sysErrorFractionBSiFlux  = 0.10; % 10%, based on a literature review
config.errorFractionPnum = 0.50; %
config.maxZeu = 200;

%% For averaging in seasons and years

config.nSeasons = 4; % fixed to 4 seasons

config.nDaysPerMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
config.cumDays = cumsum(config.nDaysPerMonth);  % [31 59 90 ... 365]

config.seasonIndices = {
    [12 1 2], ... % DJF (winter months)
    [3 4 5], ...  % MAM (spring months)
    [6 7 8], ...  % JJA (summer months)
    [9 10 11]     % SON (autumn months)
};

config.nDaysPerSeason = zeros(1, config.nSeasons);
for iSeason = 1:config.nSeasons
    config.nDaysPerSeason(iSeason) = sum(config.nDaysPerMonth(config.seasonIndices{iSeason}));
end

config.julianDayInMonth = cumsum(config.nDaysPerMonth);  

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS USED IN THIS SCRIPT
% -------------------------------------------------------------------------

% *************************************************************************

function [depths,maxNoDepths] = defineDeploymentDepths()
    
    % Define model deployment depths for sediment traps and imaging systems
    % (as in modelparameters.F90 > GetDepthsForParticleDataCollection)

    % Preallocate with a reasonably large size (adjust if needed)
    depths = [];
    currentDepth = 10;

    % 10 m intervals from 10 to 200 m
    while currentDepth <= 200
        depths(end+1) = currentDepth;
        currentDepth = currentDepth + 10;
    end

    % 50 m intervals from 250 to 1000 m
    currentDepth = 250;
    while currentDepth <= 1000
        depths(end+1) = currentDepth;
        currentDepth = currentDepth + 50;
    end

    % 100 m intervals from 1100 to 2000 m
    currentDepth = 1100;
    while currentDepth <= 2000
        depths(end+1) = currentDepth;
        currentDepth = currentDepth + 100;
    end

    % 250 m intervals from 2250 m to 5000 m
    currentDepth = 2250;
    while currentDepth <= 5000
        depths(end+1) = currentDepth;
        currentDepth = currentDepth + 250;
    end

    depths = depths'; % convert to column vector
    maxNoDepths = length(depths); % adjusted to stop at 5000 m

end % defineDeploymentDepths

% *************************************************************************

end % loadModelConfigurationParameters