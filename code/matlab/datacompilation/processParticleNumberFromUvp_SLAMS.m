% ======================================================================= %
%                                                                         %
% This script reads in particle concentration data from the UVP5          %
% instrument downloaded from the EcoTaxa repository. The data are in      %
% units of # part L-1.                                                    %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 22 May 2024                                   %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./data/raw/'));
addpath(genpath('./data/processed/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Sampled depths declarations
ECOTAXA_VERTICAL_STEP = 5; % 5 m
uvpDepths = (2.5:ECOTAXA_VERTICAL_STEP:2005)'; % m
nObsDepths = numel(uvpDepths);

% Particle size parameter declarations
NUM_SIZE_CLASSES = 15;
diameterClassEdges = zeros(NUM_SIZE_CLASSES,1);
diameterClassEdges(1) = 1; % 1 um
for i = 2:NUM_SIZE_CLASSES
    diameterClassEdges(i) = 2*diameterClassEdges(i-1);
end

% Enter the coordinates that we have used to define our locations in the
% EcoTaxa's website map
NUM_LOCS = 6;
LAT_UPPER = zeros(NUM_LOCS,1);
LAT_LOWER = zeros(NUM_LOCS,1);
LON_RIGHT = zeros(NUM_LOCS,1);
LON_LEFT = zeros(NUM_LOCS,1);

% EqPac                % OSP                   % PAP-SO               
LAT_UPPER(1) = 4;      LAT_UPPER(2) = 51.5;    LAT_UPPER(3) = 49.5;   
LAT_LOWER(1) = -4;     LAT_LOWER(2) = 49.5;    LAT_LOWER(3) = 48.5;   
LON_RIGHT(1) = -148;   LON_RIGHT(2) = -144;    LON_RIGHT(3) = -16;    
LON_LEFT(1) = -152;    LON_LEFT(2) = -146;     LON_LEFT(3) = -17;     

% BATS/OFP             % HOT/ALOHA             % HAUSGARTEN  
LAT_UPPER(4) = 32;     LAT_UPPER(5) = 23;      LAT_UPPER(6) = 80;
LAT_LOWER(4) = 29.5;   LAT_LOWER(5) = 22;      LAT_LOWER(6) = 78;
LON_RIGHT(4) = -62;    LON_RIGHT(5) = -157.5;  LON_RIGHT(6) = 5.5;
LON_LEFT(4) = -65;     LON_LEFT(5) = -158.5;   LON_LEFT(6) = 3.5;

% EcoTaxa folder definitions
SUFFIX_ECOTAXA_FOLDER_NAME = {'EqPac','OSP','PAPSO','BATSOFP','HOTALOHA','HAUSGARTEN'};

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - READ IN ECOTAXA PARTICLE FILES 
% -------------------------------------------------------------------------

% Read in the *PART* files downloaded from the EcoTaxa website
readEcoTaxaParticleFiles(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS,uvpDepths);

% Extract particle number data and classify by month
castMonthlyDistrib = calculateNumberOfCasts(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS);
maxNumCastsPerMonth = max(castMonthlyDistrib,[],'all');
uvpPnumByCast = NaN(maxNumCastsPerMonth,NUM_SIZE_CLASSES,nObsDepths,12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
 
    % The EcoTaxa data set
    load(fullfile('.','data','raw','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},'_particle_concentration.mat')),...
        'ET','sizeClassLabels')

    % Convert datetime columns and add 'month' and 'year' columns
    ET.date = datetime(ET.yyyy_mm_ddHh_mm,'format','yyyy-MM-dd');
    ET.dateString = datestr(ET.date, 'yyyy-mm-dd HH:MM:SS');
    ET.month = month(ET.date);
    ET.year = year(ET.date);

    % Classify particle data by month, depth and cast
    for iMonth = 1:12
        monthFilter = ET.month == iMonth;

        if (sum(monthFilter) > 0)
            for iDepth = 1:nObsDepths
                depthFilter = ET.Depth_m_ == uvpDepths(iDepth) & monthFilter;

                if (sum(depthFilter) > 0)
                    
                    % Crop particle number
                    pNumArray = zeros(height(ET),NUM_SIZE_CLASSES);
                    for i = 1:NUM_SIZE_CLASSES
                        pNumArray(:,i) = ET.(sizeClassLabels{i}); % # part. L-1
                    end
                    pNumArray(isnan(pNumArray)) = 0;
                    
                    % Store the filtered data in the output array
                    nInstances = sum(depthFilter);
                    uvpPnumByCast((1:nInstances),:,iDepth,iMonth,iLoc) = pNumArray(depthFilter,:);
                   
                end
                
            end % iDepth
        end
    end % iMonth
end % iLoc

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - SAVE THE DATA 
% -------------------------------------------------------------------------

save(fullfile('.','data','processed','pnum_16sc_compilation_slams.mat'),...
    'uvpPnumByCast','diameterClassEdges','uvpDepths','castMonthlyDistrib')

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function readEcoTaxaParticleFiles(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS,targetDepths)
    
for iLoc = 1:NUM_LOCS
 
    listLocalEcoTaxaDirs = getLocalEcoTaxaDirectoryNames(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS);
    pathDir = fullfile('.','data','raw','UVP5',listLocalEcoTaxaDirs{iLoc});

    % Metadata
    metadataFile = dir(fullfile(pathDir,'*metadata*.tsv'));
    M = readtable(fullfile(pathDir,metadataFile.name),...
        'FileType','text','Delimiter','tab','TreatAsEmpty',{'N/A','n/a'}); % metadata table
    
    % Keep only relevant metadata columns
    M = M(:, {'Profile', 'Cruise', 'Longitude', 'Latitude'});

    listOfParticleFiles = dir(fullfile(pathDir,'*PAR*.tsv')); % 'PAR' stands for 'particle'
    nSitesSampledForParticles = length(listOfParticleFiles);

    for iSite = 1:nSitesSampledForParticles

         thisParticleFile = listOfParticleFiles(iSite).name;
         T = readtable(fullfile(pathDir,thisParticleFile),...
             'FileType','text','Delimiter','tab','TreatAsEmpty',{'N/A','n/a'}); % file table

         % On the first iteration, identify particle size classes
         if (iSite == 1)
            idxsSizeClasses = strncmp(T.Properties.VariableNames,'LPM_',length('LPM_')); % 'LPM' stands for particle size class
            sizeClassLabels = T.Properties.VariableNames(idxsSizeClasses)'; 
         end

         % Crop relevant particle information
         Tt3 = T(:,idxsSizeClasses); % crop file table
         Tt2 = T(:,(3:5)); % date, depth and sampled volume
         Tt1 = repmat(M(iSite,:),[height(Tt3) 1]); % repeat profile, cruise, longitude and latitude data

         % Combine metadata, particle info, and particle size data
         P = [Tt1,Tt2,Tt3]; 

         if (iSite == 1 && ~isempty(P))
            vars = P.Properties.VariableNames;
            pnumData = cell2table(cell(0,length(vars)), 'VariableNames', vars);
            pnumData = P;
         elseif (iSite > 1 && ~isempty(P))
            % Append data from the current file to the existing particle data
            pnumData = [pnumData; P];
         end

    end % iLocSite
    
    ET = pnumData;

    % Extract data only for the depths of interest
    matchingDepthRows = ismember(ET.Depth_m_, targetDepths);
    ET = ET(matchingDepthRows,:);
    
    save(fullfile('.','data','raw','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},'_particle_concentration.mat')),...
        'ET','sizeClassLabels','vars')

end % iLoc

end % readEcoTaxaParticleFiles

% *************************************************************************

function listLocalEcoTaxaDirs = getLocalEcoTaxaDirectoryNames(...
    SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS)

% EcoTaxa directories are named using a specific structure consisting of a
% "prefix", "date of data download" and "suffix". This section manages
% data downloaded on various dates and organises the list of EcoTaxa
% directories according to the order specified by the parameter
% SUFFIX_ECOTAXA_FOLDER_NAME.

% Get a list of all directories in the search path and filter out non-directory entries
pathEcoTaxaDirs = dir(fullfile('.','data','raw','UVP5'));
pathEcoTaxaDirs = pathEcoTaxaDirs([pathEcoTaxaDirs.isdir]);
nameEcoTaxaDirs = {pathEcoTaxaDirs(3:end).name}; % remove '.' and '..' from the list of directories

% Create a mapping from each suffix to its index in the suffix array
suffixToIndexMap = containers.Map(SUFFIX_ECOTAXA_FOLDER_NAME, 1:NUM_LOCS);

% Initialise an array to store the reordered directory names
listLocalEcoTaxaDirs = cell(NUM_LOCS,1);
for i = 1:numel(nameEcoTaxaDirs)
    splitName = strsplit(nameEcoTaxaDirs{i}, '_');
    thisSubdirSuffix = splitName{4}; % get fourth part
    
    % Find the index of the suffix in the suffix array using the mapping
    desiredIndex = suffixToIndexMap(thisSubdirSuffix);
    
    % Place the subdirectory name in the desired position in the reordered array
    listLocalEcoTaxaDirs{desiredIndex} = nameEcoTaxaDirs{i};
end

end % getLocalEcoTaxaDirectoryNames

% *************************************************************************

function castMonthlyDistrib = calculateNumberOfCasts(SUFFIX_ECOTAXA_FOLDER_NAME,NUM_LOCS)

castMonthlyDistrib = zeros(12,NUM_LOCS);

for iLoc = 1:NUM_LOCS
 
    % The EcoTaxa data set
    load(fullfile('.','data','raw','UVP5',...
        strcat(SUFFIX_ECOTAXA_FOLDER_NAME{iLoc},'_particle_concentration.mat')),'ET')
    
    % Convert datetime columns and add 'month' and 'year' columns
    ET.date = datetime(ET.yyyy_mm_ddHh_mm,'format','yyyy-MM-dd');
    ET.dateString = datestr(ET.date, 'yyyy-mm-dd HH:MM:SS');

    for iMonth = 1:12
        monthFilter = month(ET.date) == iMonth;
        if (sum(monthFilter) > 0)
            % Get unique combinations of latitude, longitude and time
            [uniqueCombinations, ~, ~] = unique(ET{monthFilter,... 
                {'Latitude', 'Longitude','dateString'}}, 'rows');
            castMonthlyDistrib(iMonth,iLoc) = size(uniqueCombinations, 1);
        end
    end
    
end % iLoc

end % calculateNumberOfCasts

% *************************************************************************