function selectGridDomain(fullpathInputDataDir, fullpathOutputDataDir,...
    fullpathConfigDir, filenameInputTemplateGrid, filenameOutputRunGrid,...
    listLocalLats, listLocalLons, choiceTypeGridDomain, choiceZoomFactor)

% SELECTGRIDDOMAIN Creates the grid file that will be used to run SLAMS.
%
%   Written by A. Rufas, University of Oxford
%   Version 1.0 - Completed 17 Jan 2024                                  
%
% Folder paths:
%   fullpathInputDataDir  = './data/raw/';
%   fullpathOutputDataDir = './tests/test_global_esa/modelinputdata/';
%   fullpathConfigDir     = './config/';
%
% Input/output files:
%   filenameInputTemplateGrid = 'grid_MITgcm_2p8deg.mat'; 
%   filenameOutputRunGrid     = 'grid_run.mat';
%
% Choices:
%   listLocalLats = []
%   listLocalLons = []
%   choiceTypeGridDomain = 1; 
%       1: global
%       2: local
%       3: global with zoom in factor
%   choiceZoomFactor = 1;
%       any number: magnifying factor between coordinates
% 
% addpath(fullpathInputDataDir);
% addpath(fullpathOutputDataDir);
% addpath(fullpathConfigDir);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Output file
fullpathOutputFile = fullfile(fullpathOutputDataDir,filenameOutputRunGrid);

% Load template grid
load(strcat(fullpathInputDataDir,filenameInputTemplateGrid),...
    'ixBb','iyBb','Xbb','Ybb','Zsb','Zsb3d','x','y')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - MODIFY TEMPLATE GRID FILE ACCORDING TO CHOICE AND SAVE
% -------------------------------------------------------------------------

switch choiceTypeGridDomain

    case 1 % Global: use template grid as it is

        nLocs = length(Ybb);
        matchDepths = NaN(height(Zsb),nLocs);
        for iLoc = 1:nLocs
            matchDepths(:,iLoc) = Zsb3d(ixBb(iLoc),iyBb(iLoc),:);
        end
        
        clear Zsb % it includes locations that are land
        Zsb = matchDepths; % only sea locations
        save(fullpathOutputFile,'Xbb','Ybb','Zsb','ixBb','iyBb','x','y');

    case 2 % Local: extract specified locations

        fileLonsID = fopen(fullfile(fullpathConfigDir,listLocalLons),'r');
        locLons = fscanf(fileLonsID,'%f\n');
        fclose(fileLonsID); 
        
        fileLatsID = fopen(fullfile(fullpathConfigDir,listLocalLats),'r');
        locLats = fscanf(fileLatsID,'%f\n');
        fclose(fileLatsID);
        
        % Initialise arrays for specific locations
        nLocs = length(locLons);
        matchLons = NaN(nLocs,1);
        matchLats = NaN(nLocs,1);
        matchDepths = NaN(height(Zsb),nLocs);
    
        % Loop over the specific locations
        for i = 1:nLocs
            % Find the closest grid point in the global data
            [~, lonIdx] = min(abs(Xbb - locLons(i)));
            [~, latIdx] = min(abs(Ybb - locLats(i)));

            % Extract the data for the specific location
            matchLons(i) = Xbb(lonIdx);
            matchLats(i) = Ybb(latIdx);
            matchDepths(:,i) = Zsb3d(ixBb(lonIdx),iyBb(latIdx),:);
        end
        
        clear ixBb iyBb Xbb Ybb Zsb Zsb3d x y
        Xbb = matchLons;
        Ybb = matchLats;
        Zsb = matchDepths;
        save(fullpathOutputFile,'Xbb','Ybb','Zsb');
    
    case 3 % Global with zoom in factor
        
        [nr,nc] = size(Zsb);
        coordMagnification = (x(2)-x(1))*choiceZoomFactor;
        
        [lonsTemplate,latsTemplate] = ndgrid(x,y);
        [lonsZoom,latsZoom] = ndgrid(min(x):coordMagnification:max(x),...
                                     min(y):coordMagnification:max(y)); 

        seafloorDepthsTemplate = zeros(length(x),length(y));
        for iLon = 1:length(x)    
            for iLat = 1:length(y)
                if (~isnan(Zsb3d(iLon,iLat,1)))
                    iSfValidDepths = find(~isnan(Zsb3d(iLon,iLat,:)), 1, 'last');
                    seafloorDepthsTemplate(iLon,iLat) = Zsb3d(iLon,iLat,iSfValidDepths);
                end
            end
        end
        
        seafloorDepthZoom = interpn(lonsTemplate,latsTemplate,seafloorDepthsTemplate,lonsZoom,latsZoom,'linear');

        %         figure(1)
        %         pcolor(flipud(rot90(seafloorDepthsTemplate)))
        %         shading interp
        %         colormap(jet)
        %         box on
        %         colorbar
        %         title('Seafloor depths template grid')
        %
        %         figure(2)
        %         pcolor(flipud(rot90(seafloorDepthZoom)))
        %         shading interp
        %         colormap(jet)
        %         box on
        %         colorbar
        %         title('Seafloor depths zoom grid')

        xZoom = lonsZoom(:,1);
        yZoom = latsZoom(1,:)';
        nxZoom = length(xZoom);
        nyZoom = length(yZoom);
        Zsb3dZoom = NaN(nxZoom,nyZoom,nr); % NaN as default value for unused depths
        ixBbZoom = zeros(nxZoom*nyZoom,1); 
        iyBbZoom = zeros(nxZoom*nyZoom,1);
        XbbZoom = zeros(nxZoom*nyZoom,1);
        YbbZoom = zeros(nxZoom*nyZoom,1);
        ZsbZoom = NaN(nr,nxZoom*nyZoom,1); % NaN as default value for unused depths

        % Track the number of locations with non-empty depths
        iLoc = 0;
        for iLat = 1:nyZoom
            for iLon = 1:nxZoom
                % Compute the depths for the current location
                depthsZoom = (5:10:seafloorDepthZoom(iLon,iLat));
                
                if ~isempty(depthsZoom)
                    % Increment location counter
                    iLoc = iLoc + 1;
                    
                    % Store depths in the 3D matrix and 2D matrix
                    Zsb3dZoom(iLon,iLat,(1:length(depthsZoom))) = depthsZoom;
                    ZsbZoom((1:length(depthsZoom)),iLoc) = depthsZoom;
                    
                    % Store the indices and coordinates
                    ixBbZoom(iLoc) = iLon;
                    iyBbZoom(iLoc) = iLat;
                    XbbZoom(iLoc) = xZoom(iLon);
                    YbbZoom(iLoc) = yZoom(iLat);
                end        
            end
        end

        % Truncate arrays to only include completed locations
        nCompletedLocs = iLoc; % number of completed locations (non-empty depths)
        ixBbZoom = ixBbZoom(1:nCompletedLocs);
        iyBbZoom = iyBbZoom(1:nCompletedLocs);
        XbbZoom = XbbZoom(1:nCompletedLocs);
        YbbZoom = YbbZoom(1:nCompletedLocs);
        ZsbZoom = ZsbZoom(:,1:nCompletedLocs);

        % Save the outputs
        clear ixBb iyBb Xbb Ybb Zsb Zsb3d x y
        ixBb = ixBbZoom;
        iyBb = iyBbZoom;
        Xbb = XbbZoom;
        Ybb = YbbZoom;
        Zsb = ZsbZoom;
        Zsb3d = Zsb3dZoom;
        x = xZoom;
        y = yZoom;

        save(fullpathOutputFile,'x','y','ixBb','iyBb','Xbb','Ybb','Zsb3d','Zsb');

end

end % selectGridDomain
