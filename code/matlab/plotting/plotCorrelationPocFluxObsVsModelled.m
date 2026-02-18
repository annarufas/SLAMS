function plotCorrelationPocFluxObsVsModelled(figureSubfolderName,...
    figurePrefixName,modFluxMonthlyDh,obsFluxDataMonthly,bicepExFluxMonthly,...
    filenameTimeseriesInformation)

% Load station-related information from observation compilation
load(fullfile('.','data','interim',filenameTimeseriesInformation),'STATION_NAMES')

% Rearrange model locations to match observational order
desiredLocationOrder = {'HOT/ALOHA','BATS/OFP','EqPac','PAP-SO','OSP','HAUSGARTEN'};
currentModLocationOrder = STATION_NAMES;
[~,reorderModLocIdx] = ismember(desiredLocationOrder,currentModLocationOrder); % get reordering indices
nLocs = length(STATION_NAMES);

% Extract POC flux data
iTracer = 1;
modPocFlux    = squeeze(modFluxMonthlyDh(1,:,reorderModLocIdx,:,iTracer)); % 1 = export depth
modPocFlux    = mean(modPocFlux,3,'omitnan'); % 12 x 6
obsPocFlux    = squeeze(obsFluxDataMonthly{iTracer}(1,:,reorderModLocIdx)); % 12 x 6
bicepExFluxMonthly_perm = permute(bicepExFluxMonthly, [2 1 3]); % 6 x 12 --> 12 x 6
bicepLiPocFlux          = squeeze(bicepExFluxMonthly_perm(:,reorderModLocIdx,3)); % 6 x 12
bicepHensonPocFlux      = squeeze(bicepExFluxMonthly_perm(:,reorderModLocIdx,2)); 
bicepDunnePocFlux       = squeeze(bicepExFluxMonthly_perm(:,reorderModLocIdx,1)); 

% Reshape (12 x 6 --> 72 x 1)
modPocFlux_vec = modPocFlux(:); % 72 x 1 column vector
obsPocFlux_vec = obsPocFlux(:);      
bicepLiPocFlux_vec = bicepLiPocFlux(:);
bicepHensonPocFlux_vec = bicepHensonPocFlux(:);
bicepDunnePocFlux_vec = bicepDunnePocFlux(:);

% Matching location IDs: repeat each location ID for all months
locIDs = repelem((1:nLocs)',12); % 72 x 1

%% First plot: observations vs SLAMS-2.0

% Measurement sets for looping
methods    = {obsPocFlux_vec, bicepDunnePocFlux_vec, bicepHensonPocFlux_vec, bicepLiPocFlux_vec};
methodNames= {'{\it In situ}','Satellite (BICEP-Dunne)','Satellite (BICEP-Henson)','Satellite (BICEP-Li)'};

% Colour map for locations
coloursLocs = brewermap(nLocs,'RdYlBu');

% Determine global axis limits across all comparisons
allVals = [modPocFlux_vec; obsPocFlux_vec; bicepDunnePocFlux_vec; bicepHensonPocFlux_vec; bicepLiPocFlux_vec];
lims = [min(allVals) max(allVals)];

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.40],'Color','w')
haxis = zeros(2,2);

for iType = 1:4
    haxis(iType) = subaxis(2,2,iType,'Spacing',0.07,'Padding',0,'Margin',0.09);
    ax(iType).pos = get(haxis(iType),'Position');
    
    % Shift all plots up
    if (iType == 2 || iType == 4)
        ax(iType).pos(1) = ax(iType).pos(1)-0.075;
    else
        ax(iType).pos(1) = ax(iType).pos(1)-0.025; 
    end
    if (iType == 1 || iType == 2)
        ax(iType).pos(2) = ax(iType).pos(2)+0.030;
    else
        ax(iType).pos(2) = ax(iType).pos(2)+0.005;
    end
    set(haxis(iType),'Position',ax(iType).pos)
    
    Y = modPocFlux_vec;
    X = methods{iType};
    locVec = locIDs(:); % 

    % Find indices where neither X nor Y is NaN
    validIdx = ~isnan(X) & ~isnan(Y);
    X = X(validIdx);
    Y = Y(validIdx);
    locVec = locVec(validIdx); % filter locations as well
  
    % Scatter points coloured by location
    hold on;
    for iLoc = 1:nLocs
        idx = locVec == iLoc;
        scatter(haxis(iType), X(idx), Y(idx), 60, 'o', 'MarkerEdgeColor','k',...
            'MarkerFaceColor', coloursLocs(iLoc,:), 'LineWidth',0.5);
    end

    % Axis equal and square
    axis square;
    xlim([0 280]); % xlim(lims);
    ylim([0 440]);

    % Set tick labels every 100 units
    xLimits = xlim; 
    yLimits = ylim;

    xticksValues = xLimits(1):100:xLimits(2);
    xticks(xticksValues);

    yticksValues = yLimits(1):100:yLimits(2);
    yticks(yticksValues);    
    
    % 1:1 line
    plot(haxis(iType), lims, lims, 'k-', 'LineWidth', 1);
    
    % Regression line
    regressCoeff = polyfit(X, Y, 1);
    slope = regressCoeff(1);
    xfit = linspace(lims(1), lims(2), 100);
    yfit = polyval(regressCoeff, xfit);
    plot(haxis(iType), xfit, yfit, 'k--', 'LineWidth', 1.5);
    
    % Performance statistics: get correlation coefficient, bias (me) and RMSE
    [corrCoeff, pvalue] = corr(X,Y,'type','Spearman');
    [rmse,~] = calcRootMeanSquaredError(X,Y); % greek letter: phi
    [me,~] = calcMeanError(X,Y); % aka bias, greek letter: delta
    N = length(X);

    str = sprintf('\\rho = %.2f\nslope = %.2f\n\\delta = %.0f\nRMSE = %.0f', ...
        corrCoeff,slope,me,rmse);

    if pvalue <= 0.05
        % Position near upper right corner
        if iType ~= 1
            xt = max(xlim) - 0.02*max(xlim); 
            yt = max(ylim) - 0.02*max(ylim);
            text(xt,yt,str,'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top');
        % Position text near bottom right corner
        else
            xt = max(xlim) - 0.02*range(xlim);
            yt = min(ylim) + 0.02*range(ylim);
            text(xt,yt,str,'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom');
        end
    end

    title(methodNames{iType},'FontSize',13)

end

% Create an invisible axes spanning the entire figure
han = axes(gcf, 'Visible', 'off');

% Turn on visibility just for the labels
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

% Shared labels
xlabel(han, 'Observation of POC export flux (mg C m^{-2} d^{-1})', 'FontSize', 12);
ylabel(han, 'SLAMS-2.0 estimate of POC export flux (mg C m^{-2} d^{-1})', 'FontSize', 12);

xl = han.XLabel;
xl.Position(2) = xl.Position(2) - 0.03;  % negative moves down, positive moves up

yl = han.YLabel;
yl.Position(1) = yl.Position(1) - 0.07;  % negative moves left, positive moves right

% Legend
lg = legend(haxis(1),desiredLocationOrder);
lg.Position(1) = 0.83;
lg.Position(2) = 0.765;
lg.ItemTokenSize = [11,4];
lg.FontSize = 11;
lg.NumColumns = 1;
set(lg, 'Box', 'off');

% Save
saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_','obs_vs_model'))

%% Second plot: in situ vs satellite-based estimates

% Measurement sets for looping
methods    = {bicepDunnePocFlux_vec, bicepHensonPocFlux_vec, bicepLiPocFlux_vec};
methodNames= {'BICEP-Dunne','BICEP-Henson','BICEP-Li'};

% Colour map for locations
coloursLocs = brewermap(nLocs,'RdYlBu');

% Determine global axis limits across all comparisons
allVals = [obsPocFlux_vec; bicepLiPocFlux_vec; bicepHensonPocFlux_vec; bicepDunnePocFlux_vec];
lims = [min(allVals) max(allVals)];

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.50 0.23],'Color','w')
haxis = zeros(1,3);

for iType = 1:3
    haxis(iType) = subaxis(1,3,iType,'Spacing',0.04,'Padding',0,'Margin',0.14);
    ax(iType).pos = get(haxis(iType),'Position');

    % Shift all plots to the left and up
    if (iType == 1)
        ax(iType).pos(1) = ax(iType).pos(1)-0.060;
    elseif (iType == 2)
        ax(iType).pos(1) = ax(iType).pos(1)-0.050;
    elseif (iType == 3)
        ax(iType).pos(1) = ax(iType).pos(1)-0.040;
    end
    ax(iType).pos(2) = ax(iType).pos(2)+0.030;
    set(haxis(iType),'Position',ax(iType).pos)

    X = obsPocFlux_vec;
    Y = methods{iType};
    locVec = locIDs(:); 

    % Find indices where neither X nor Y is NaN
    validIdx = ~isnan(X) & ~isnan(Y);
    X = X(validIdx);
    Y = Y(validIdx);
    locVec = locVec(validIdx); % filter locations as well

    % Scatter points coloured by location
    hold on;
    for iLoc = 1:nLocs
        idx = locVec == iLoc;
        scatter(haxis(iType), X(idx), Y(idx), 60, 'o', 'MarkerEdgeColor','k',...
            'MarkerFaceColor', coloursLocs(iLoc,:), 'LineWidth',0.5);
    end

    % Axis equal and square
    axis square;
    xlim([0 300]); % xlim(lims);
    ylim([0 200]);

    % Set tick labels every 100 units
    xLimits = xlim;  % get current limits
    yLimits = ylim;

    xticksValues = xLimits(1):100:xLimits(2);
    xticks(xticksValues);

    yticksValues = yLimits(1):100:yLimits(2);
    yticks(yticksValues);    

    % 1:1 line
    plot(haxis(iType), lims, lims, 'k-', 'LineWidth', 1);
    
    % Regression line
    regressCoeff = polyfit(X, Y, 1);
    slope = regressCoeff(1);
    xfit = linspace(lims(1), lims(2), 100);
    yfit = polyval(regressCoeff, xfit);
    plot(haxis(iType), xfit, yfit, 'k--', 'LineWidth', 1.5);
    
    % Performance statistics: get correlation coefficient, bias (me) and RMSE
    [corrCoeff, pvalue] = corr(X,Y,'type','Spearman');
    [rmse,~] = calcRootMeanSquaredError(X,Y); % greek letter: phi
    [me,~] = calcMeanError(X,Y); % aka bias, greek letter: delta
    N = length(X);

    % Position near upper right corner
    xt = max(xlim)-0.02*max(xlim); 
    yt = max(ylim)-0.02*max(ylim);

    if (pvalue <= 0.05)
        str = sprintf('\\rho = %.2f\nslope = %.2f\n\\delta = %.0f\nRMSE = %.0f', ...
                  corrCoeff, slope, me, rmse);
        text(xt,yt,str,'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','top');
    end

    if iType == 2
        xlabel('{\it In situ} observation of POC export flux (mg C m^{-2} d^{-1})','FontSize',12);
    end
    if iType == 1
        ylabel('Satellite-based estimate (mg C m^{-2} d^{-1})','FontSize',12);
    end
    title(methodNames{iType},'FontSize',13)

end

% Legend
lg = legend(desiredLocationOrder);
lg.Position(1) = 0.84;
lg.Position(2) = 0.59;
lg.ItemTokenSize = [11,4];
lg.FontSize = 11;
lg.NumColumns = 1;
set(lg, 'Box', 'off');

% Save
saveFigureInFolder(figureSubfolderName,strcat(figurePrefixName,'_','insitu_vs_satell'))

end % plotCorrelationPocFluxObsVsModelled
