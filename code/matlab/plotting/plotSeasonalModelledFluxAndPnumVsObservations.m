function plotSeasonalModelledFluxAndPnumVsObservations(figureSubfolderName,...
    figurePrefixName,config,modFlux,obsFlux,obsError,bicepExFlux,bicepExDepth,...
    modPnum,uvpPnum,uvpDepths,filenameTimeseriesInformation,filenameObsPocFlux,...
    filenameObsPicAndBsiFlux,filenameObsPnumUvp5)
%%
% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')

% Rearrange locations to match desired order
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
currentLocationOrder = STATION_NAMES;
[~,reorderLocIdx] = ismember(desiredLocationOrder,currentLocationOrder); % get reordering indices
nLocs = length(STATION_NAMES);

% POC flux compilation observed depths
load(fullfile('.','data','processed',filenameObsPocFlux),'LOC_DEPTH_HORIZONS_POC_OBS_ADAPT')

% PIC and bSi flux compilation observed depths
load(fullfile('.','data','processed',filenameObsPicAndBsiFlux),'LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT')

seasonLabel = {'Winter','Spring','Summer','Autumn'};
idxSeason = config.seasonIndices;

bicepModelTags = {'Dunne','Henson','Li'};
myColourSchemeBicep = brewermap(numel(bicepModelTags),'*YlOrRd');

%% Plot biogeochemical fluxes

for iTracer = 1:3

    modProfileValues = modFlux(:,:,:,iTracer);
    obsDhValues = obsFlux{iTracer}(:,:,:);
    obsDhErr = obsError{iTracer}(:,:,:);

    if (iTracer == 1)
        depthsTmp = LOC_DEPTH_HORIZONS_POC_OBS_ADAPT; % size: [12, 6, 2, 4]
    elseif (iTracer == 2)
        depthsTmp = squeeze(LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(:,:,:,:,1)); 
    elseif (iTracer == 3)
        depthsTmp = squeeze(LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(:,:,:,:,2)); 
    end

    depthsTmpSeasonal = NaN(4,nLocs,2,4);
    for iSeason = 1:4
        seasonalMean = mean(depthsTmp(idxSeason{iSeason},:,:,:),1,'omitnan');  % size: [4, 6, 2, 4]
        depthsTmpSeasonal(iSeason,:,:,:) = squeeze(seasonalMean);
    end
    obsDhDepths = squeeze(mean(depthsTmpSeasonal(:,:,:,:),3)); % mean of upper and lower

    modSedTrapProfileDepths = config.availSedTrapDeployDepths;
    obsDhDepths = permute(obsDhDepths, [3 1 2]); % reorder dimensions

    % Extract modelled data
    modProfileValuesForPlotting = NaN(size(modProfileValues));
    modProfileDepthsForPlotting = NaN(size(modProfileValues));
    for iSeason = 1:4
        for iLoc = 1:nLocs
            modProfileValuesForPlotting(:,iSeason,iLoc) = modProfileValues(:,iSeason,iLoc);
            modProfileDepthsForPlotting(:,iSeason,iLoc) = modSedTrapProfileDepths(:);    
        end
    end

    % Apply the right index location order
    obsDhValues = obsDhValues(:,:,reorderLocIdx);
    obsDhErr =  obsDhErr(:,:,reorderLocIdx);
    obsDhDepths = obsDhDepths(:,:,reorderLocIdx);
    modProfileValuesForPlotting = modProfileValuesForPlotting(:,:,reorderLocIdx);
    modProfileDepthsForPlotting = modProfileDepthsForPlotting(:,:,reorderLocIdx);
    bicepExFlux  = bicepExFlux(reorderLocIdx,:,:);
    bicepExDepth = bicepExDepth(:,reorderLocIdx);

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.50],'Color','w')
    haxis = zeros(nLocs*4,1);

    for iSeason = 1:4
        for iLoc = 1:nLocs

            % Compute linear subplot index (row-wise fill)
            iSubplot = (iSeason - 1) * 6 + iLoc;
            haxis(iSubplot) = subaxis(4,nLocs,iSubplot,'Spacing',0.038,'Padding',0,'Margin',0.11);
            ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    
            %%%%%%%%%%%%%% MODELLED DATA %%%%%%%%%%%%%%
            
            currSeasonTotPnumModData = squeeze(modProfileValuesForPlotting(:,iSeason,iLoc)); % mg m-2 d-1
            modSedTrapProfileDepths = squeeze(modProfileDepthsForPlotting(:,iSeason,iLoc));
    
            %%%%%%%%%%%%%% OBSERVED DATA FROM COMPILATION %%%%%%%%%%%%%%
    
            currSeasonTotPnumObsData = squeeze(obsDhValues(:,iSeason,iLoc)); % mg m-2 d-1
            currSeasonObsDataErr     = squeeze(obsDhErr(:,iSeason,iLoc)); 
            currSeasonObsDepths      = squeeze(obsDhDepths(:,iSeason,iLoc)); 

            %%%%%%%%%%%%%% OBSERVED DATA FROM BICEP %%%%%%%%%%%%%%

            if (iTracer == 1)
                currSeasonBicepData  = squeeze(bicepExFlux(iLoc,iSeason,:)); % mg m-2 d-1
                currSeasonBicepDepth = squeeze(bicepExDepth(iSeason,iLoc));
            end

            % .................................................................

            % Modelled fluxes
            hpl(1) = plot(currSeasonTotPnumModData, modSedTrapProfileDepths, 'o',...
                'MarkerEdgeColor', [0.698, 0.875, 0.541], 'MarkerFaceColor', [0.698, 0.875, 0.541],...
                'LineStyle', 'none', 'LineWidth', 0.5,...
                'DisplayName', 'Modelled');
            hold on
    
            % Observations from compilation
            hpl(2) = plot(currSeasonTotPnumObsData, currSeasonObsDepths, 'o',...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.106, 0.471, 0.216],...
                'LineStyle', 'none', 'LineWidth', 0.5,...
                'DisplayName', 'Observed (trap & radionuclide)');
            hold on

            % Add horizontal error bars manually
            for i = 1:length(currSeasonTotPnumObsData)
                line([currSeasonTotPnumObsData(i)-currSeasonObsDataErr(i), currSeasonTotPnumObsData(i)+currSeasonObsDataErr(i)], [currSeasonObsDepths(i), currSeasonObsDepths(i)], 'Color', 'k')
            end

            % Observations from BICEP export flux
            if (iTracer == 1) 
                for iModel = 1:3
                     hpl(iModel+2) = plot(currSeasonBicepData(iModel), currSeasonBicepDepth, 'o',...
                        'MarkerEdgeColor', 'k', 'MarkerFaceColor', myColourSchemeBicep(iModel,:),...
                        'LineStyle', 'none', 'LineWidth', 0.5,...
                        'DisplayName', strcat('Observed (BICEP-',bicepModelTags{iModel},')'));
                    hold on
                end
            end

            hold on
            box on

            % Fine tune for some locations
            if iTracer == 1
                if iLoc == 1 % HOT/ALOHA
                    xlim([0 105])
                    xTickValues = 0:50:100;
                elseif iLoc == 2 % BATS/OFP 
                    xlim([0 175])
                    xTickValues = 0:75:150;
                elseif iLoc == 3 % EqPac
                    xlim([0 380])
                    xTickValues = 0:150:300;
                elseif iLoc == 4 % PAP-SO
                    xlim([0 440])
                    xTickValues = 0:200:400;       
                elseif iLoc == 5 % OSP
                    xlim([0 280])
                    xTickValues = 0:100:200;
                elseif iLoc == 6
                    xlim([0 110])
                    xTickValues = 0:50:100;
                end
            elseif iTracer == 2 
                if iLoc == 1 % HOT/ALOHA
                    xlim([0 180])
                    xTickValues = 0:75:150;
                elseif iLoc == 2 % BATS/OFP
                    xlim([0 180])
                    xTickValues = 0:75:150;
                elseif iLoc == 3 % EqPac
                    xlim([0 420])
                    xTickValues = 0:200:400;
                elseif iLoc == 4 % PAP-SO
                    xlim([0 550])
                    xTickValues = 0:250:500;
                elseif iLoc == 5 % OSP
                    xlim([0 420])
                    xTickValues = 0:200:400;
                elseif iLoc == 6 % HAUSGARTEN
                    xlim([0 130])
                    xTickValues = 0:50:100;
                end
            elseif iTracer == 3 
                if iLoc == 1 % HOT/ALOHA
                    xlim([0 10])
                    xTickValues = 0:5:10;
                elseif iLoc == 2
                    xlim([0 36])
                    xTickValues = 0:15:30;
                elseif iLoc == 3
                    xlim([0 52])
                    xTickValues = 0:25:50;
                elseif iLoc == 4
                    xlim([0 300])
                    xTickValues = 0:100:200;
                elseif iLoc == 5
                    xlim([0 170])
                    xTickValues = 0:75:150;
                elseif iLoc == 6
                    xlim([0 35])
                    xTickValues = 0:15:30;
                end
            end

            xticks(xTickValues)
            xticklabels(xTickValues)

            ylim([0 1600])
            yticks([0 500 1000 1500])
            % set(gca,'Yscale','log')
            axh = gca;
            axh.YAxis.TickDirection = 'out';
            axh.TickLength = [0.03, 0.03]; % make tick marks longer
    
            if (iLoc == 1)
                yticklabels({'0','500','1000','1500'})
            else
                yticklabels([])
            end
    
            set(gca,'YDir','Reverse','XAxisLocation','Bottom','xlabel',[],'ylabel',[],'FontSize', 10)
    
            if iSeason == 1
                text(0.5, 1.35, desiredLocationOrder{iLoc}, ...
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
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

        end % iLoc
    end % iSeason

    % Add season labels above each row
    for iSeason = 1:4
        % Get subplot indices for the current row
        rowIdx = (iSeason - 1) * nLocs + (1:nLocs);
    
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

        % % Some touches to the legend
        % if (iTracer == 1)
        %     lg = legend([hpl(1) hpl(2) hpl(3) hpl(4) hpl(5)]);
        %     lg.Position(1) = 0.37; lg.Position(2) = 0.005;
        %     lg.Orientation = 'horizontal';
        %     lg.ItemTokenSize = [20,50];
        %     lg.FontSize = 11;
        %     lg.NumColumns = 1;
        %     lg.Box = 'on';
        % end
        % 

    % Give common xlabel, ylabel and title to your figure
    % Create a new axis
    a = axes;
    if (iTracer == 1)
        xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',14);
    elseif (iTracer == 2)
        xl = xlabel('PIC flux (mg calcite m^{-2} d^{-1})','FontSize',14);
    elseif (iTracer == 3)
        xl = xlabel('bSi flux (mg opal m^{-2} d^{-1})','FontSize',14);
    end
    yl = ylabel('Depth (m)','FontSize',14);
    
    % Specify visibility of the current axis as 'off'
    a.Visible = 'off';
    % Specify visibility of Title, XLabel, and YLabel as 'on'
    xl.Visible = 'on';
    yl.Visible = 'on';
    yl.Position(1) = yl.Position(1) - 0.07; 
    xl.Position(1) = 0.5; xl.Position(2) = xl.Position(2) - 0.045;

    if (iTracer == 1)
        plotTag = 'pocflux';
    elseif (iTracer == 2)
        plotTag = 'picflux';
    elseif (iTracer == 3)
        plotTag = 'bsiflux';
    end
    saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_',plotTag))

end % iTracer

%% Plot particle numbers

% Apply the right index location order
uvpPnum = uvpPnum(:,:,:,reorderLocIdx);
uvpDepths = uvpDepths(:,:,:,reorderLocIdx);
modPnum = modPnum(:,:,:,reorderLocIdx);
modImagSysProfileDepths = config.availImagSysDeployDepths;

% Modelled size classes
modPartSizeClassEdges = config.particleDiameterClassBounds;

% Obserevd size classes
load(fullfile('.','data','processed',filenameObsPnumUvp5),'diameterClassEdges')
obsPartSizeClassEdges = diameterClassEdges;

% Indices in modPartSizeClassEdges that match obs
[isMatch, idxMod] = ismember(obsPartSizeClassEdges, modPartSizeClassEdges);
modIdx = idxMod(isMatch);

% Colour palettes
coloursMod = flipud(brewermap(3,'PuBuGn'));
coloursObs = flipud(brewermap(3,'YlOrRd'));

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.50],'Color','w')
haxis = zeros(nLocs*4,1);
hLegend = gobjects(6,1);

for iSeason = 1:4
    for iLoc = 1:nLocs

        % Compute linear subplot index (row-wise fill)
        iSubplot = (iSeason - 1) * nLocs + iLoc;
        haxis(iSubplot) = subaxis(4,nLocs,iSubplot,'Spacing',0.038,'Padding',0,'Margin',0.11);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    
        %%%%%%%%%%%%%% MODELLED DATA %%%%%%%%%%%%%%
            
        currSeasonTotPnumModData = squeeze(sum(modPnum(modIdx,:,iSeason,iLoc),1,'omitnan'));
        currSeasonSmallPnumModData = squeeze(sum(modPnum(2:7,:,iSeason,iLoc),1,'omitnan'));
        currSeasonLargePnumModData = squeeze(sum(modPnum(8:end,:,iSeason,iLoc),1,'omitnan'));

        % Discard last modelled depth (seafloor)
        idxsSampledDepths = find(currSeasonTotPnumModData ~= 0);
        xtot = currSeasonTotPnumModData(1:numel(idxsSampledDepths)-1)';
        xsmall = currSeasonSmallPnumModData(1:numel(idxsSampledDepths)-1)';
        xlarge = currSeasonLargePnumModData(1:numel(idxsSampledDepths)-1)';
        y = modImagSysProfileDepths(1:numel(idxsSampledDepths)-1);
        
        %%%%%%%%%%%%%% OBSERVED DATA %%%%%%%%%%%%%%

        currSeasonTotPnumObsData   = squeeze(sum(uvpPnum(:,:,iSeason,iLoc),1,'omitnan'))'; 
        currSeasonSmallPnumObsData   = squeeze(sum(uvpPnum(1:6,:,iSeason,iLoc),1,'omitnan'))'; 
        currSeasonLargePnumObsData   = squeeze(sum(uvpPnum(7:end,:,iSeason,iLoc),1,'omitnan'))'; 
        currSeasonObsDepths = squeeze(uvpDepths(1,:,iSeason,iLoc))'; 

        % .................................................................

        isLegendAxis = (iSeason == 1 && iLoc == 1);

        % Modelled pnum
        h2 = plot(xsmall, y, 'o',...
            'MarkerEdgeColor', coloursMod(2,:), 'MarkerFaceColor', coloursMod(1,:),...
            'LineStyle', 'none', 'LineWidth', 0.5); hold on;
        h3 = plot(xlarge, y, 'o',...
            'MarkerEdgeColor', coloursMod(3,:), 'MarkerFaceColor', coloursMod(3,:),...
            'LineStyle', 'none', 'LineWidth', 0.5); hold on;

        % Observed pnum (UVP5)
        h5 = plot(currSeasonSmallPnumObsData, currSeasonObsDepths, 'o',...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', coloursObs(1,:),...
            'LineStyle', 'none', 'LineWidth', 0.5); hold on;
        h6 = plot(currSeasonLargePnumObsData, currSeasonObsDepths, 'o',...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', coloursObs(3,:),...
            'LineStyle', 'none', 'LineWidth', 0.5); hold off;

        box on

        if isLegendAxis
            hLegend = [h2 h3 h5 h6];
        end

        set(gca, 'XScale', 'log')
        xlim([1 5e7])                                      
        xticks([1 1e3 1e6])                         
        xticklabels({'10^0','10^3','10^6'})   
        set(gca, 'XTickLabelRotation', 0) % prevent rotation

        ylim([0 1600])
        yticks([0 500 1000 1500])
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.03, 0.03]; % make tick marks longer

        if (iLoc == 1)
            yticklabels({'0','500','1000','1500'})
        else
            yticklabels([])
        end

        set(gca,'YDir','Reverse','XAxisLocation','Bottom','xlabel',[],'ylabel',[],'FontSize', 11)

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
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1)-0.030; 
        set(haxis(iSubplot),'Position',ax(iSubplot).pos)

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
    rowIdx = (iSeason - 1) * 6 + (1:6);

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

% Some touches to the legend
lg = legend(hLegend, ...
    {'Model small','Model large', ...
     'UVP5 small','UVP5 large'});
lg.Position(1) = 0.88;
lg.Position(2) = 0.07;
lg.Orientation = 'vertical';
lg.ItemTokenSize = [20,50];
lg.FontSize = 11;
lg.NumColumns = 1;
lg.Box = 'off';

% Give common xlabel, ylabel and title to your figure
% Create a new axis
a = axes;
xl = xlabel('Particle number concentration (# L^{-1})','FontSize',14);
yl = ylabel('Depth (m)','FontSize',14);

% Specify visibility of the current axis as 'off'
a.Visible = 'off';
% Specify visibility of Title, XLabel, and YLabel as 'on'
xl.Visible = 'on';
yl.Visible = 'on';
yl.Position(1) = yl.Position(1) - 0.1; 
xl.Position(1) = 0.45; xl.Position(2) = xl.Position(2) - 0.05;
    
saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_','pnum'))

end % plotSeasonalModelledFluxAndPnumVsObservations