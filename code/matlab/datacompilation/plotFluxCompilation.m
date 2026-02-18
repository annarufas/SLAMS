
% ======================================================================= %
%                                                                         %
% This script plots the POC, PIC and bSi flux compilation of sediment     %
% trap and radionuclide data. It has been modified from the script        %
% "plotPocFluxFromTrapAndRadCompilation.m".                               %
%                                                                         %
% The script has 4 sections:                                              %
%   Section 1 - Presets.                                                  %
%   Section 2 - Plot POC, PIC and bSi flux data by station and month.     %
%   Section 3 - Plot no. entries by month, location and depth horizon     %
%               for POC, PIC and bSi flux).                               %
%   Section 4 - Plot monthly fluxes and their error by location and depth %
%               horizon.                                                  %
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed 25 Apr 2025                                   %
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

filenameInputPocFluxCompilation       = 'pocflux_compilation_slams.mat';
filenameInputPicAndBsiFluxCompilation = 'picandbsiflux_compilation_slams.mat';
filenameInputTimeseriesInformation    = 'timeseries_station_information_slams.mat';

load(fullfile('.','data','processed',filenameInputPocFluxCompilation),...
    'pocFluxRawProfileValues','pocFluxRawProfileDepths','pocFluxRawProfileDataType',...
    'pocFluxRawDhValues_cell','pocFluxRawDhDepths_cell','pocFluxRawDhTag_cell','pocFluxRawDhDataType_cell',...
    'pocFluxMonthlyDhAvg','pocFluxMonthlyDhN','pocFluxMonthlyDhErrTot')

load(fullfile('.','data','processed',filenameInputPicAndBsiFluxCompilation),...
    'picbsiFluxRawProfileValues','picbsiFluxRawProfileDepths','picbsiFluxRawProfileDataType',...
    'picbsiFluxRawDhValues_cell','picbsiFluxRawDhDepths_cell','picbsiFluxRawDhTag_cell','picbsiFluxRawDhDataType_cell',...
    'picbsiFluxMonthlyDhAvg','picbsiFluxMonthlyDhN','picbsiFluxMonthlyDhErrTot')

load(fullfile('.','data','interim',filenameInputTimeseriesInformation),...
    'LOC_DEPTH_HORIZONS','STATION_NAMES','STATION_TAGS','MAX_NUM_VALUES_PER_MONTH')

% Parameters
NUM_TARGET_DEPTHS = size(LOC_DEPTH_HORIZONS,4);
NUM_LOCS = length(STATION_NAMES);
MOLAR_MASS_CARBON = 12.011; % g mol-1
MOLAR_MASS_SI = 28.0855; % g mol-1
monthLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
myColours = parula(NUM_TARGET_DEPTHS);

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 2 - PLOT FLUX DATA BY STATION AND MONTH
% -------------------------------------------------------------------------

for iTracer = 1:3

    if (iTracer == 1)
        rawProfileValues = pocFluxRawProfileValues;
        rawProfileDepths = pocFluxRawProfileDepths;
        rawProfileDataType = pocFluxRawProfileDataType;
        rawDhValues_cell = pocFluxRawDhValues_cell;
        rawDhDepths_cell = pocFluxRawDhDepths_cell;
        rawDhDataType_cell = pocFluxRawDhDataType_cell;
        rawDhTag_cell = pocFluxRawDhTag_cell;
    elseif (iTracer == 2)
        rawProfileValues = picbsiFluxRawProfileValues(:,:,:,1);
        rawProfileDepths = picbsiFluxRawProfileDepths(:,:,:,1);
        rawProfileDataType = picbsiFluxRawProfileDataType(:,:,:,1);
        rawDhValues_cell = picbsiFluxRawDhValues_cell(:,:,:,1);
        rawDhDepths_cell = picbsiFluxRawDhDepths_cell(:,:,:,1);
        rawDhDataType_cell = picbsiFluxRawDhDataType_cell(:,:,:,1);
        rawDhTag_cell = picbsiFluxRawDhTag_cell(:,:,:,1);
    elseif (iTracer == 3)
        rawProfileValues = picbsiFluxRawProfileValues(:,:,:,2);
        rawProfileDepths = picbsiFluxRawProfileDepths(:,:,:,2);
        rawProfileDataType = picbsiFluxRawProfileDataType(:,:,:,2);
        rawDhValues_cell = picbsiFluxRawDhValues_cell(:,:,:,2);
        rawDhDepths_cell = picbsiFluxRawDhDepths_cell(:,:,:,2);
        rawDhDataType_cell = picbsiFluxRawDhDataType_cell(:,:,:,2);
        rawDhTag_cell = picbsiFluxRawDhTag_cell(:,:,:,2);
    end

    for iLoc = 1:NUM_LOCS
    
        figure()
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.75],'Color','w')
        haxis = zeros(12,1);
        
        for iMonth = 1:12
    
            haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.055,'Padding',0,'Margin',0.13);
            ax(iMonth).pos = get(haxis(iMonth),'Position');
    
            %%%%%%%%%%%%%% ALL DATA %%%%%%%%%%%%%%
            
            if (iTracer == 3)
                yall = MOLAR_MASS_SI.*squeeze(rawProfileValues(:,iMonth,iLoc)); % mmol m-2 d-1 --> mg m-2 d-1
            else
                yall = MOLAR_MASS_CARBON.*squeeze(rawProfileValues(:,iMonth,iLoc)); % mmol m-2 d-1 --> mg m-2 d-1
            end
            xall = squeeze(rawProfileDepths(:,iMonth,iLoc));
            tyall = squeeze(rawProfileDataType(:,iMonth,iLoc)); % trap vs rad
    
            %%%%%%%%%%%%%% DH DATA %%%%%%%%%%%%%%
    
            currMonthDhData     = zeros(NUM_TARGET_DEPTHS,MAX_NUM_VALUES_PER_MONTH);
            currMonthDhDataType = cell(NUM_TARGET_DEPTHS,MAX_NUM_VALUES_PER_MONTH);
            currMonthDhTag      = cell(NUM_TARGET_DEPTHS,MAX_NUM_VALUES_PER_MONTH);
            currMonthDhDepths   = zeros(NUM_TARGET_DEPTHS,MAX_NUM_VALUES_PER_MONTH);
            currMonthDhN        = zeros(NUM_TARGET_DEPTHS,1);
    
            for iDh = 1:NUM_TARGET_DEPTHS
                allmyvals      = rawDhValues_cell{iDh,iMonth,iLoc};
                allmydepths    = rawDhDepths_cell{iDh,iMonth,iLoc};
                allmydatatypes = rawDhDataType_cell{iDh,iMonth,iLoc};
                allmydhtags    = rawDhTag_cell{iDh,iMonth,iLoc};
                if ~isempty(allmyvals)
                    vals = str2num(allmyvals);
                    depths = str2num(allmydepths);
                    types = textscan(allmydatatypes,'%s');
                    tags = textscan(allmydhtags,'%s');
                    % Store processed data
                    currMonthDhN(iDh) = numel(vals);
                    if (iTracer == 3)
                        currMonthDhData(iDh,1:currMonthDhN(iDh)) = MOLAR_MASS_SI.*vals; % mmol m-2 d-1 --> mg m-2 d-1
                    else
                        currMonthDhData(iDh,1:currMonthDhN(iDh)) = MOLAR_MASS_CARBON.*vals; % mmol m-2 d-1 --> mg m-2 d-1
                    end
                    currMonthDhDepths(iDh,1:currMonthDhN(iDh)) = depths;
                    % Store data types and tags
                    for iDataPoint = 1:numel(vals)
                        currMonthDhDataType{iDh,iDataPoint} = types{1}{iDataPoint};
                        currMonthDhTag{iDh,iDataPoint} = tags{1}{iDataPoint};
                    end
                end
            end % iDh
    
            [rowIdxs,colIdxs,ydh] = find(currMonthDhData);
            nDataPointsInDh = numel(ydh);
            xdh = zeros(nDataPointsInDh,1);
            tydh = cell(nDataPointsInDh,1);
            for iDataPoint = 1:nDataPointsInDh
                xdh(iDataPoint) = currMonthDhDepths(rowIdxs(iDataPoint),colIdxs(iDataPoint));
                tydh(iDataPoint) = currMonthDhTag(rowIdxs(iDataPoint),colIdxs(iDataPoint));
            end
            
            % .................................................................
    
            % We have 4 possible combinations. Mask the entries of interest.
            % If there are no entries for the specific combination 
            % (e.g., zeu-radioisotope), plot NaN, so that the plot generates 
            % an entry for that combination
                
            % All observations (grey)

            d01 = plot(yall, xall, 'o', 'MarkerEdgeColor', [0.85 0.85 0.85],...
                'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5,... 
                'HandleVisibility', 'off');
            max01 = max(yall);
            hold on
    
            % Observations by depth horizon (zeu, mesoupp, mesolow and zmeso)
    
            % zeu
            mask = strcmp(tydh,'zeu');
            if (sum(mask) == 0)
                d1 = plot(NaN, NaN, 'o','MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(4,:), 'Linewidth', 0.5,...
                    'DisplayName', 'z_{eu}');
            else
                d1 = plot(ydh(mask), xdh(mask), 'o', 'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(4,:), 'LineWidth', 0.5,...
                    'DisplayName', 'z_{eu}');
            end
            max1 = max(ydh(mask));
            hold on
            
            % mesoupp
            mask = strcmp(tydh,'mesoupp');
            if (sum(mask) == 0)
                d2 = plot(NaN, NaN, 'o','MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(3,:), 'Linewidth', 0.5,...
                    'DisplayName', 'mesoupp');
            else
                d2 = plot(ydh(mask), xdh(mask), 'o', 'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(3,:), 'LineWidth', 0.5,...
                    'DisplayName', 'mesoupp');
            end
            max2 = max(ydh(mask));
            hold on
            
            % mesolow
            mask = strcmp(tydh,'mesolow');
            if (sum(mask) == 0)
                d3 = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(2,:), 'LineWidth', 0.5,...
                    'DisplayName', 'mesolow');
            else
                d3 = plot(ydh(mask), xdh(mask), 'o', 'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(2,:), 'LineWidth', 0.5,...
                    'DisplayName', 'mesolow');
            end
            max3 = max(ydh(mask));
            hold on
            
            % zmeso
            mask = strcmp(tydh,'zmeso');
            if (sum(mask) == 0)
                d4 = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(1,:), 'LineWidth', 0.5,...
                    'DisplayName', 'z_{meso}');
            else
                d4 = plot(ydh(mask), xdh(mask), 'o', 'MarkerEdgeColor', 'k',...
                    'MarkerFaceColor', myColours(1,:), 'LineWidth', 0.5,...
                    'DisplayName', 'z_{meso}');
            end
            max4 = max(ydh(mask));
            hold off
    
            box on
            
            if (iTracer == 1)
                if (iLoc == 1) % EqPac
                    xlim([0 500])
                    xTickValues = 0:200:400;
                elseif (iLoc == 2) % OSP
                    xlim([0 200])
                    xTickValues = 0:75:150;       
                elseif (iLoc == 3) % PAP-SO
                    xlim([0 380])
                    xTickValues = 0:150:300;
                elseif (iLoc == 4) % BATS/OFP
                    xlim([0 350])
                    xTickValues = 0:150:300;
                elseif (iLoc == 5) % HOT/ALOHA
                    xlim([0 100])
                    xTickValues = 0:50:100;
                elseif (iLoc == 6) % HAUSGARTEN
                    xlim([0 70])
                    xTickValues = 0:25:50;
                end
            elseif (iTracer == 2)
                if (iLoc == 1) % EqPac
                    xlim([0 20])
                    xTickValues = 0:10:20;
                elseif (iLoc == 2) % OSP
                    xlim([0 78])
                    xTickValues = 0:30:60;       
                elseif (iLoc == 3) % PAP-SO
                    xlim([0 50])
                    xTickValues = 0:25:50;
                elseif (iLoc == 4) % BATS/OFP
                    xlim([0 22])
                    xTickValues = 0:10:20;
                elseif (iLoc == 5) % HOT/ALOHA
                    xlim([0 11])
                    xTickValues = 0:5:10;
                elseif (iLoc == 6) % HAUSGARTEN
                    xlim([0 16])
                    xTickValues = 0:5:15;
                end
            elseif (iTracer == 3)
                if (iLoc == 1) % EqPac
                    xlim([0 24])
                    xTickValues = 0:10:20;
                elseif (iLoc == 2) % OSP
                    xlim([0 150])
                    xTickValues = 0:75:150;       
                elseif (iLoc == 3) % PAP-SO
                    xlim([0 280])
                    xTickValues = 0:100:200;
                elseif (iLoc == 4) % BATS/OFP
                    xlim([0 15])
                    xTickValues = 0:7:14;
                elseif (iLoc == 5) % HOT/ALOHA
                    xlim([0 10])
                    xTickValues = 0:5:10;
                elseif (iLoc == 6) % HAUSGARTEN
                    xlim([0 21])
                    xTickValues = 0:10:20;
                end
            end

            xticks(xTickValues)
            xticklabels(xTickValues)
    
            ylim([15 2000])
            yticks([20 100 500 1000 2000])
            set(gca,'Yscale','log')
            axh = gca;
            axh.YAxis.TickDirection = 'out';
            axh.TickLength = [0.03, 0.03]; % make tick marks longer
    
            if (iMonth == 1 || iMonth == 4 || iMonth ==7 || iMonth == 10)
                yticklabels({'20','100','500','1000','2000'})
            else
                yticklabels([])
            end
    
            set(gca,'YDir','Reverse','XAxisLocation','Bottom','xlabel',[],'ylabel',[],'FontSize', 11)
    
            title(monthLabel(iMonth),'FontSize',16)
    
        end % iMonth
    
        % Shift all plots a little bit to the right and up
        ax(1).pos(1) = ax(1).pos(1)+0.040; ax(1).pos(2) = ax(1).pos(2)+0.040; 
        ax(2).pos(1) = ax(2).pos(1)+0.025; ax(2).pos(2) = ax(2).pos(2)+0.040; 
        ax(3).pos(1) = ax(3).pos(1)+0.010; ax(3).pos(2) = ax(3).pos(2)+0.040;  
        ax(4).pos(1) = ax(4).pos(1)+0.040; ax(4).pos(2) = ax(4).pos(2)+0.040; 
        ax(5).pos(1) = ax(5).pos(1)+0.025; ax(5).pos(2) = ax(5).pos(2)+0.040; 
        ax(6).pos(1) = ax(6).pos(1)+0.010; ax(6).pos(2) = ax(6).pos(2)+0.040; 
        ax(7).pos(1) = ax(7).pos(1)+0.040; ax(7).pos(2) = ax(7).pos(2)+0.040; 
        ax(8).pos(1) = ax(8).pos(1)+0.025; ax(8).pos(2) = ax(8).pos(2)+0.040; 
        ax(9).pos(1) = ax(9).pos(1)+0.010; ax(9).pos(2) = ax(9).pos(2)+0.040; 
        ax(10).pos(1) = ax(10).pos(1)+0.040; ax(10).pos(2) = ax(10).pos(2)+0.040; 
        ax(11).pos(1) = ax(11).pos(1)+0.025; ax(11).pos(2) = ax(11).pos(2)+0.040; 
        ax(12).pos(1) = ax(12).pos(1)+0.010; ax(12).pos(2) = ax(12).pos(2)+0.040; 
        for iMonth = 1:12
            set(haxis(iMonth),'Position',ax(iMonth).pos) 
        end
    
        % Some touches to the legend
        lg = legend([d1 d2 d3 d4]);
        lg.Position(1) = 0.35; lg.Position(2) = -0.005;
        lg.Orientation = 'horizontal';
        lg.ItemTokenSize = [20,50];
        lg.FontSize = 11;
        lg.NumColumns = 2;
        lg.Box = 'on';
    
        % Give common xlabel, ylabel and title to your figure
        % Create a new axis
        a = axes;
        t = title(STATION_NAMES(iLoc),'FontSize',18);
        if (iTracer == 1)
            xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',18);
        elseif (iTracer == 2)
            xl = xlabel('PIC flux (mg C m^{-2} d^{-1})','FontSize',18);
        elseif (iTracer == 3)
            xl = xlabel('bSi flux (mg Si m^{-2} d^{-1})','FontSize',18);
        end
        
        yl = ylabel('Depth (m)','FontSize',18);
        % Specify visibility of the current axis as 'off'
        a.Visible = 'off';
        % Specify visibility of Title, XLabel, and YLabel as 'on'
        t.Visible = 'on';
        xl.Visible = 'on';
        yl.Visible = 'on';
        yl.Position(1) = yl.Position(1) - 0.005; yl.Position(2) = yl.Position(2) + 0.02; 
        xl.Position(1) = 0.5; xl.Position(2) = xl.Position(2) + 0.065;
        t.Position(1) = t.Position(1); t.Position(2) = t.Position(2) + 0.035;
        
        if (iTracer == 1)
            plotTag = 'pocflux';
        elseif (iTracer == 2)
            plotTag = 'picflux';
        elseif (iTracer == 3)
            plotTag = 'bsiflux';
        end
        saveFigureInFolder('compilation',strcat('compilation_',plotTag,'_',STATION_TAGS{iLoc}))
    
    end % iLoc
end % iFluxTracer

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 3 - PLOT FIGURE SHOWING NO. ENTRIES BY MONTH, LOCATION AND DEPTH
% HORIZON
% -------------------------------------------------------------------------

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.60],'Color','w')
haxis = zeros(3,1);

for iTracer = 1:3

    if (iTracer == 1)
        theNumberOfDataPoints = pocFluxMonthlyDhN(:,:,:); % nDepths x 12 months x NUM_LOCS
    elseif (iTracer == 2)
        theNumberOfDataPoints = picbsiFluxMonthlyDhN(:,:,:,1); 
    elseif (iTracer == 3)
        theNumberOfDataPoints = picbsiFluxMonthlyDhN(:,:,:,2); 
    end

    theNumberOfDataPoints_permutted = permute(theNumberOfDataPoints, [2 3 1]); % 12 months x NUM_LOCS x nDepths
    
    % Swap locations: currently, the order is (1) EqPac, (2) OSP, (3) PAP-SO, 
    % (4) BATS/OFP, (5) HOT/ALOHA and (6) HAUSGARTEN, and the desired order is 
    % (1) HOT/ALOHA, (2) BATS/OFP, (3) EqPac, (4) PAP-SO, (5) OSP and (6) HAUSGARTEN.
    theNumberOfDataPoints_swapped = theNumberOfDataPoints_permutted;
    theNumberOfDataPoints_swapped(:, [3, 5, 4, 2, 1, 6], :) = theNumberOfDataPoints_permutted(:, [1, 2, 3, 4, 5, 6], :);
    
    % Number of data points by depth
    summedData = sum(sum(theNumberOfDataPoints, 2), 3); 
    if (iTracer == 1)
        fprintf('\nPOC flux:')
    elseif (iTracer == 2)
        fprintf('\nPIC flux:')
    elseif (iTracer == 3)
        fprintf('\nbSi flux:')
    end
    fprintf('\n%0.1f%% data points at zeu', 100.*(summedData(1)/sum(summedData)))
    fprintf('\n%0.1f%% data points at mesoupp', 100.*(summedData(2)/sum(summedData)))
    fprintf('\n%0.1f%% data points at mesolow', 100.*(summedData(3)/sum(summedData)))
    fprintf('\n%0.1f%% data points at zmeso', 100.*(summedData(4)/sum(summedData)))
    
    % .........................................................................
    
    haxis(iTracer) = subaxis(3,1,iTracer,'Spacing',0.08,'Padding',0,'Margin',0.10);
    ax(iTracer).pos = get(haxis(iTracer),'Position');
    
    yy = flipdim(theNumberOfDataPoints_swapped,3); % to have zmeso at the bottom of the plot instead of at the top
    h = plotBarStackGroups(yy, monthLabel); % plot groups of stacked bars

    % Change the colors of each bar segment
    coloursBarSegments = parula(size(h,2));
    coloursBarSegments = repelem(coloursBarSegments,size(h,1),1); 
    coloursBarSegments = mat2cell(coloursBarSegments,ones(size(coloursBarSegments,1),1),3);
    set(h,{'FaceColor'},coloursBarSegments)
    
    if (iTracer == 1)
        ylim([0 200])
        yTickValues = 0:50:200;
    elseif (iTracer == 2)
        ylim([0 100])
        yTickValues = 0:25:100;
    elseif (iTracer == 3)
        ylim([0 50])
        yTickValues = 0:10:50;
    end
    yticks(yTickValues)
    yticklabels(yTickValues)

    yl = ylabel('Number of data points');
    yl.Position(1) = yl.Position(1) - 0.5;
    box on

    if (iTracer == 1)
        title('POC flux')
    elseif (iTracer == 2)
        title('PIC flux');
    elseif (iTracer == 3)
        title('bSi flux');
    end

    % Grid lines
    ax = gca;
    ax.YGrid = 'on'; % horizontal grid lines
    ax.XGrid = 'off'; % no vertical grid lines
    set(gca, 'TickLength', [0 0])
    
    % Legend
    lg = legend('Base of mesopelagic zone','Lower mesopelagic','Upper mesopelagic','Base of euphotic zone');
    lg.Position(1) = 0.38; lg.Position(2) = -0.020;
    lg.Orientation = 'horizontal';
    set(lg,'Box','off') 
    
    set(gca, 'FontSize', 14);

    clear ax

end % iTracer

saveFigureInFolder('compilation','compilation_numberdatapoints')

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 4 - PLOT MONTHLY FLUXES AND THEIR ERROR BY LOCATION AND DEPTH 
% HORIZON
% -------------------------------------------------------------------------

flipMyColours = flipud(myColours); % flip upside down

% Transform units of observed data, as used in previous plots:
% mmol C m-2 d-1 --> mg C m-2 d-1 and
% mmol Si m-2 d-1 --> mg Si m-2 d-1 
obsData = {MOLAR_MASS_CARBON.*pocFluxMonthlyDhAvg,...
           MOLAR_MASS_CARBON.*picbsiFluxMonthlyDhAvg(:,:,:,1),...
           MOLAR_MASS_SI.*picbsiFluxMonthlyDhAvg(:,:,:,2)};
obsError = {MOLAR_MASS_CARBON.*pocFluxMonthlyDhErrTot,...
            MOLAR_MASS_CARBON.*picbsiFluxMonthlyDhErrTot(:,:,:,1),...
            MOLAR_MASS_SI.*picbsiFluxMonthlyDhErrTot(:,:,:,2)};

% Outlier removal in observations
[obsDataNoOutliers,obsErrorNoOutliers] = removeOutliersInFluxObservations(...
    obsData,obsError,NUM_TARGET_DEPTHS,NUM_LOCS);

for iTracer = 1:3
    valsNoOutliers = obsData{iTracer}; % obsDataNoOutliers{iTracer};
    errNoOutliers = obsError{iTracer}; % obsErrorNoOutliers{iTracer};

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.70 0.50],'Color','w') 
    haxis = zeros(NUM_TARGET_DEPTHS,NUM_LOCS);
    iSubplot = 0;
    
    for iDh = 1:NUM_TARGET_DEPTHS
        for iStation = 1:NUM_LOCS
            
            % Re-order
            switch iStation
                case 1
                    iLoc = 5; % HOT/ALOHA
                case 2
                    iLoc = 4; % BATS/OFP
                case 3
                    iLoc = 1; % EqPac
                case 4
                    iLoc = 3; % PAP-SO
                case 5
                    iLoc = 2; % OSP
                case 6
                    iLoc = 6; % HAUSGARTEN
            end
            
            iSubplot = iSubplot + 1;
    
            haxis(iSubplot) = subaxis(NUM_TARGET_DEPTHS,NUM_LOCS,iSubplot,'Spacing',0.01,'Padding',0.01,'Margin', 0.09);
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
    
            hbar = bar(haxis(iSubplot),(1:12),squeeze(valsNoOutliers(iDh,:,iLoc)),...
                'BarWidth',0.75,'FaceColor','flat');
            hbar.CData(:,:) = repmat(flipMyColours(iDh,:),[12 1]);
            hold on
            
            her = errorbar(haxis(iSubplot),(1:12),squeeze(valsNoOutliers(iDh,:,iLoc)),...
                zeros(size(squeeze(valsNoOutliers(iDh,:,iLoc)))),squeeze(errNoOutliers(iDh,:,iLoc)));    
            her.Color = [0 0 0];                            
            her.LineStyle = 'none'; 
            hold off
            
            % Base of euphotic
            if (iSubplot >= 1 && iSubplot <= 6)
                if (iTracer == 1)
                    ylim([0 300])
                    yTickValues = 0:50:300;
                elseif (iTracer == 2)
                    ylim([0 60])
                    yTickValues = 0:20:60;
                elseif (iTracer == 3)
                    ylim([0 100])
                    yTickValues = 0:25:100;
                end
            % Upper mesopelagic
            elseif (iSubplot > 6 && iSubplot <= 12)
                if (iTracer == 1)
                    ylim([0 200])
                    yTickValues = 0:50:200;
                elseif (iTracer == 2)
                    ylim([0 25])
                    yTickValues = 0:5:25;
                elseif (iTracer == 3)
                    ylim([0 100])
                    yTickValues = 0:25:100;
                end
            % Lower mesopelagic
            elseif (iSubplot > 12 && iSubplot <= 18)
                if (iTracer == 1)
                    ylim([0 100])
                    yTickValues = 0:25:100;
                elseif (iTracer == 2)
                    ylim([0 20])
                    yTickValues = 0:5:20;
                elseif (iTracer == 3)
                    ylim([0 60])
                    yTickValues = 0:20:60;
                end
            % Base of mesopelagic
            else
                if (iTracer == 1)
                    ylim([0 30])
                    yTickValues = 0:5:30;
                elseif (iTracer == 2)
                    ylim([0 15])
                    yTickValues = 0:5:15;
                elseif (iTracer == 3)
                    ylim([0 30])
                    yTickValues = 0:5:30;
                end
            end
            
            yticks(yTickValues)
            yticklabels(yTickValues)
            ytickformat('%.0f')

            if (iSubplot == 1 || iSubplot == 7 || iSubplot == 13 || iSubplot == 19)
                if (iTracer == 1)
                    yl = ylabel('POC flux (mg C m^{-2} d^{-1})');
                elseif (iTracer == 2)
                    yl = ylabel('PIC flux (mg C m^{-2} d^{-1})');
                elseif (iTracer == 3)
                    yl = ylabel('bSi flux (mg Si m^{-2} d^{-1})');
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
                tl = title(STATION_NAMES(iLoc),'FontSize',14);
                tl.Visible = 'on';
                tl.Position(2) = tl.Position(2) + 0.20;
            end
    
            grid on;
            axh.XGrid = 'off';
            axh.YGrid = 'on';
    
        end % iStation 
    end % iDh
    
    if (iTracer == 1)
        plotTag = 'pocflux';
    elseif (iTracer == 2)
        plotTag = 'picflux';
    elseif (iTracer == 3)
        plotTag = 'bsiflux';
    end
    saveFigureInFolder('compilation',strcat('compilation_flux_by_month_and_station_',plotTag))

end % iTracer