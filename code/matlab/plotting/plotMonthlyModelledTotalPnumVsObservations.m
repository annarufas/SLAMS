function plotMonthlyModelledTotalPnumVsObservations(figureSubfolderName,figurePrefixName,...
    config,modPnumMonthlyProfile,uvpPnumMonthlyTargetValues,uvpPnumMonthlyTargetDepths,...
    filenameTimeseriesInformation)

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES','STATION_TAGS')
nLocs = length(STATION_NAMES);

monthLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% Plot total particle number data 

for iLoc = 1:nLocs

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.75],'Color','w')
    haxis = zeros(12,1);
    
    for iMonth = 1:12
        haxis(iMonth) = subaxis(4,3,iMonth,'Spacing',0.055,'Padding',0,'Margin',0.13);
        ax(iMonth).pos = get(haxis(iMonth),'Position');

        %%%%%%%%%%%%%% MODELLED DATA %%%%%%%%%%%%%%
            
        currMonthModData = squeeze(sum(modPnumMonthlyProfile(:,:,iMonth,iLoc),1,'omitnan'));
        currMonthModDepths = config.availImagSysDeployDepths;

        % Discard last modelled depth (sea floor)
        idxsSampledDepths = find(currMonthModData ~= 0);
        x = currMonthModData(1:numel(idxsSampledDepths)-1)';
        y = currMonthModDepths(1:numel(idxsSampledDepths)-1);
        
        %%%%%%%%%%%%%% OBSERVED DATA %%%%%%%%%%%%%%

        currMonthObsData   = squeeze(sum(uvpPnumMonthlyTargetValues(:,:,iMonth,iLoc),1,'omitnan'))'; 
        currMonthObsDepths = squeeze(uvpPnumMonthlyTargetDepths(1,:,iMonth,iLoc))'; 

        % .................................................................

        % Modelled fluxes
        plot(x, y, 'o',...
            'MarkerEdgeColor', [0.6, 0.8, 1], 'MarkerFaceColor', [0.6, 0.8, 1], 'LineStyle', 'none', 'LineWidth', 0.5);
        hold on

        % Observations
        plot(currMonthObsData, currMonthObsDepths, 'o',...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineStyle', 'none', 'LineWidth', 0.5);
        hold off

        box on

        set(gca, 'XScale', 'log')
        xlim([1 1e10])                                      
        xticks([1 1000 1e6 1e9])                         
        xticklabels({'1','10^3','10^6','10^9'})   
        set(gca, 'XTickLabelRotation', 0)  % Prevent rotation

        ylim([0 2000])
        yticks([0 500 1000 1500 2000])
        axh = gca;
        axh.YAxis.TickDirection = 'out';
        axh.TickLength = [0.03, 0.03]; % make tick marks longer

        if (iMonth == 1 || iMonth == 4 || iMonth ==7 || iMonth == 10)
            yticklabels({'0','500','1000','1500','2000'})
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
    lg = legend('Model estimates','Observations');
    lg.Position(1) = 0.38;
    lg.Position(2) = 0.05;
    lg.Orientation = 'horizontal';
    lg.ItemTokenSize = [20,50];
    lg.FontSize = 11;
    lg.NumColumns = 2;
    lg.Box = 'on';

    % Give common xlabel, ylabel and title to your figure
    % Create a new axis
    a = axes;
    t = title(STATION_NAMES(iLoc),'FontSize',18);
    xl = xlabel('Particle number concentration (# L^{-1})','FontSize',18);

    yl = ylabel('Depth (m)','FontSize',18);
    % Specify visibility of the current axis as 'off'
    a.Visible = 'off';
    % Specify visibility of Title, XLabel, and YLabel as 'on'
    t.Visible = 'on';
    xl.Visible = 'on';
    yl.Visible = 'on';
    yl.Position(1) = yl.Position(1) - 0.005; yl.Position(2) = yl.Position(2) + 0.02; 
    xl.Position(1) = 0.5; xl.Position(2) = xl.Position(2) + 0.05;
    t.Position(1) = t.Position(1); t.Position(2) = t.Position(2) + 0.035;

    saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_',STATION_TAGS{iLoc}))

end % iLoc

end % plotMonthlyModelledTotalPnumVsObservations