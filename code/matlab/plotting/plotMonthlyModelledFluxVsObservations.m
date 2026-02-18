function plotMonthlyModelledFluxVsObservations(isNormalised,figureSubfolderName,...
    figurePrefixName,config,modFlux,obsFlux,obsError,bicepExFlux,bicepExDepth,...
    filenameTimeseriesInformation,filenameObsPocFlux,filenameObsPicAndBsiFlux)

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES','STATION_TAGS')
nLocs = length(STATION_NAMES);

% POC flux compilation observed depths
load(fullfile('.','data','processed',filenameObsPocFlux),'LOC_DEPTH_HORIZONS_POC_OBS_ADAPT')

% PIC and bSi flux compilation observed depths
load(fullfile('.','data','processed',filenameObsPicAndBsiFlux),'LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT')

% Plot
monthLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
bicepModelTags = {'Dunne','Henson','Li'};
myColourSchemeBicep = brewermap(numel(bicepModelTags),'*YlOrRd');

for iTracer = 1:3

    modProfileValues = modFlux(:,:,:,iTracer);
    obsDhValues = obsFlux{iTracer}(:,:,:);
    obsDhErr = obsError{iTracer}(:,:,:);

    if (iTracer == 1)
        obsDhDepths = squeeze(mean(LOC_DEPTH_HORIZONS_POC_OBS_ADAPT(:,:,:,:),3)); % mean of upper and lower
    elseif (iTracer == 2)
        obsDhDepths = squeeze(mean(LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(:,:,:,:,1),3)); % mean of upper and lower
    elseif (iTracer == 3)
        obsDhDepths = squeeze(mean(LOC_DEPTH_HORIZONS_PICANDBSI_OBS_ADAPT(:,:,:,:,2),3)); % mean of upper and lower
    end

    modProfileDepths = config.availSedTrapDeployDepths;
    obsDhDepths = permute(obsDhDepths, [3 1 2]); % reorder dimensions

    % Extract modelled data from the normalisation depth
    modProfileValuesForPlotting = NaN(size(modProfileValues));
    modProfileDepthsForPlotting = NaN(size(modProfileValues));

    if isNormalised

        for iMonth = 1:12
            for iLoc = 1:nLocs
    
                firstShallowestObs = obsDhValues(1,iMonth,iLoc);
                secondShallowestObs = obsDhValues(2,iMonth,iLoc);
    
                firstShallowestObsDepth = obsDhDepths(1,iMonth,iLoc);
                secondShallowestObsDepth = obsDhDepths(2,iMonth,iLoc);
    
                idxClosestModDepth = NaN;
    
                if ~isnan(firstShallowestObs)
                
                    % Find the closest depth greater than or equal to firstShallowestObsDepth
                    [~,idxClosestModDepth] = min(abs(config.availSedTrapDeployDepths - firstShallowestObsDepth));
    
                elseif (isnan(firstShallowestObs) && ~isnan(secondShallowestObs))
    
                    % Find the closest depth greater than or equal to secondShallowestObsDepth
                    [~,idxClosestModDepth] = min(abs(config.availSedTrapDeployDepths - secondShallowestObsDepth));
                end
    
                % Extract valid modelled depth range
                if ~isnan(idxClosestModDepth)
                    modProfileValuesForPlotting(idxClosestModDepth:end,iMonth,iLoc) = modProfileValues(idxClosestModDepth:end,iMonth,iLoc);
                    modProfileDepthsForPlotting(idxClosestModDepth:end,iMonth,iLoc) = modProfileDepths(idxClosestModDepth:end);    
                end
            end
        end

    else

        for iMonth = 1:12
            for iLoc = 1:nLocs
                modProfileValuesForPlotting(:,iMonth,iLoc) = modProfileValues(:,iMonth,iLoc);
                modProfileDepthsForPlotting(:,iMonth,iLoc) = modProfileDepths(:);    
            end
        end

    end % isNormalised

    for iLoc = 1:nLocs
    
        figure()
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.75],'Color','w')
        haxis = zeros(12,1);
        
        for iMonth = 1:12
    
            haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.055,'Padding',0,'Margin',0.13);
            ax(iMonth).pos = get(haxis(iMonth),'Position');
    
            %%%%%%%%%%%%%% MODELLED DATA %%%%%%%%%%%%%%
            
            currMonthModData = squeeze(modProfileValuesForPlotting(:,iMonth,iLoc)); % mg m-2 d-1
            currMonthModDepths = squeeze(modProfileDepthsForPlotting(:,iMonth,iLoc));
    
            %%%%%%%%%%%%%% OBSERVED DATA FROM COMPILATION %%%%%%%%%%%%%%
    
            currMonthObsData    = squeeze(obsDhValues(:,iMonth,iLoc)); % mg m-2 d-1
            currMonthObsDataErr = squeeze(obsDhErr(:,iMonth,iLoc)); 
            currMonthObsDepths  = squeeze(obsDhDepths(:,iMonth,iLoc)); 

            %%%%%%%%%%%%%% OBSERVED DATA FROM BICEP %%%%%%%%%%%%%%

            if (iTracer == 1)
                currMonthBicepData  = squeeze(bicepExFlux(iLoc,iMonth,:)); % mg m-2 d-1
                currMonthBicepDepth = squeeze(bicepExDepth(iMonth,iLoc));
            end

            % .................................................................

            if isNormalised

                % Only plot those profiles if there are >= 2 observations 
                if sum(~isnan(currMonthObsData)) >= 2
    
                    % Modelled fluxes
                    plot(currMonthModData, currMonthModDepths, 'o',...
                        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineStyle', 'none', 'LineWidth', 0.5);
                    hold on
            
                    % Observations
                    plot(currMonthObsData, currMonthObsDepths, 'o',...
                        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineStyle', 'none', 'LineWidth', 0.5);
                
                end

            else

                % Modelled fluxes
                hpl(1) = plot(currMonthModData, currMonthModDepths, 'o',...
                    'MarkerEdgeColor', [0 0.7 0], 'MarkerFaceColor', [0 0.7 0],...
                    'LineStyle', 'none', 'LineWidth', 0.5,...
                    'DisplayName', 'Modelled');
                hold on
        
                % Observations from compilation
                hpl(2) = plot(currMonthObsData, currMonthObsDepths, 'o',...
                    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k',...
                    'LineStyle', 'none', 'LineWidth', 0.5,...
                    'DisplayName', 'Observed (trap & radionuclide)');
                hold on

                % Add horizontal error bars manually
                for i = 1:length(currMonthObsData)
                    line([currMonthObsData(i)-currMonthObsDataErr(i), currMonthObsData(i)+currMonthObsDataErr(i)], [currMonthObsDepths(i), currMonthObsDepths(i)], 'Color', 'k')
                end

                if (iTracer == 1)
                    % Observations from BICEP export flux
                    for iModel = 1:3
                         hpl(iModel+2) = plot(currMonthBicepData(iModel), currMonthBicepDepth, 'o',...
                            'MarkerEdgeColor', 'k', 'MarkerFaceColor', myColourSchemeBicep(iModel,:),...
                            'LineStyle', 'none', 'LineWidth', 0.5,...
                            'DisplayName', strcat('Observed (BICEP-',bicepModelTags{iModel},')'));
                        hold on
                    end
                end

            end
            hold on

            box on

            % Deal with axis labelling    
            if ~isNormalised

                if (iTracer == 1)
                    xlim([0 160])
                    xTickValues = 0:75:150;
                elseif (iTracer == 2)
                    xlim([0 400])
                    xTickValues = 0:150:300;
                elseif (iTracer == 3)
                    xlim([0 70])
                    xTickValues = 0:30:60;
                end

                % Fine tune for some locations
                if iTracer == 1
                    if iLoc == 1 % EqPac
                        xlim([0 400])
                        xTickValues = 0:150:300;
                    elseif iLoc == 2 % OSP
                        xlim([0 270])
                        xTickValues = 0:100:200;       
                    elseif iLoc == 3 % PAP-SO
                        xlim([0 520])
                        xTickValues = 0:200:400;
                    elseif iLoc == 4 % BATS
                        xlim([0 210])
                        xTickValues = 0:100:200;  
                    end
                elseif iTracer == 2  
                    if iLoc == 4 || iLoc == 5
                        xlim([0 280])
                        xTickValues = 0:100:200;
                    elseif iLoc == 2 
                        xlim([0 500])
                        xTickValues = 0:200:400;
                    elseif iLoc == 3 
                        xlim([0 700])
                        xTickValues = 0:300:600;
                    elseif iLoc == 6
                        xlim([0 160])
                        xTickValues = 0:75:150;
                    end
                elseif iTracer == 3 
                    if iLoc == 1 % EqPac
                        xlim([0 90])
                        xTickValues = 0:40:80;
                    elseif iLoc == 2 % OSP
                        xlim([0 180])
                        xTickValues = 0:75:150;
                    elseif iLoc == 3
                        xlim([0 140])
                        xTickValues = 0:60:120;
                    elseif iLoc == 4
                        xlim([0 38])
                        xTickValues = 0:15:30;
                    elseif iLoc == 5
                        xlim([0 10])
                        xTickValues = 0:5:10;
                    end
                end

            else

                if (iTracer == 1)
                    xlim([0 1.6])
                    xTickValues = 0:0.5:1.5;
                elseif (iTracer == 2)
                    xlim([0 3.5])
                    xTickValues = 0:1:3;
                elseif (iTracer == 3)
                    xlim([0 1.6])
                    xTickValues = 0:0.5:1.5;
                end

                % Fine tune for some locations
                if (iTracer == 1 && iLoc == 6)
                    xlim([0 30])
                    xTickValues = 0:10:30;
                elseif (iTracer == 1 && iLoc == 2)
                    xlim([0 3])
                    xTickValues = 0:1:3;
                elseif (iTracer == 2 && (iLoc == 3 || iLoc == 4))
                    xlim([0 9])
                    xTickValues = 0:4:8;
                elseif (iTracer == 2 && iLoc == 6)
                    xlim([0 250])
                    xTickValues = 0:100:200;
                elseif (iTracer == 3 && iLoc == 3)
                    xlim([0 8])
                    xTickValues = 0:3:6;    
                elseif (iTracer == 3 && (iLoc == 2 || iLoc == 4))
                    xlim([0 14])
                    xTickValues = 0:5:10;
                elseif (iTracer == 3 && iLoc == 6)
                    xlim([0 220])
                    xTickValues = 0:100:200;
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
    
            if (iMonth == 1 || iMonth == 4 || iMonth ==7 || iMonth == 10)
                yticklabels({'0','500','1000','1500'})
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
        if (iTracer == 1)
            lg = legend([hpl(1) hpl(2) hpl(3) hpl(4) hpl(5)]);
            lg.Position(1) = 0.37; lg.Position(2) = 0.005;
            lg.Orientation = 'horizontal';
            lg.ItemTokenSize = [20,50];
            lg.FontSize = 11;
            lg.NumColumns = 1;
            lg.Box = 'on';
        end
    
        % Give common xlabel, ylabel and title to your figure
        % Create a new axis
        a = axes;
        t = title(STATION_NAMES(iLoc),'FontSize',18);

        if ~isNormalised
            if (iTracer == 1)
                xl = xlabel('POC flux (mg C m^{-2} d^{-1})','FontSize',18);
            elseif (iTracer == 2)
                xl = xlabel('PIC flux (mg calcite m^{-2} d^{-1})','FontSize',18);
            elseif (iTracer == 3)
                xl = xlabel('bSi flux (mg opal m^{-2} d^{-1})','FontSize',18);
            end
        else
            if (iTracer == 1)
                xl = xlabel('Normalised POC flux','FontSize',18);
            elseif (iTracer == 2)
                xl = xlabel('Normalised PIC flux','FontSize',18);
            elseif (iTracer == 3)
                xl = xlabel('Normalised bSi flux','FontSize',18);
            end
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
        saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_',plotTag,'_',STATION_TAGS{iLoc}))

    end % iLoc
end % iTracer

end % plotMonthlyModelledFluxVsObservations