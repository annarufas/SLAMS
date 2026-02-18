function [metrics,labelMetricNames] = calculatePerformanceMetrics(insituData,modelData)

    % Fit a linear model
    regressCoeff = polyfit(insituData, modelData, 1);
    slope = regressCoeff(1);

    % Correlation coefficient
    [corrCoeff,pvalue] = corr(insituData, modelData, 'Type', 'Spearman');
    rSquared = corrCoeff^2;
    
    % Error metrics
    [rmse,~] = calcRootMeanSquaredError(insituData,modelData); % greek letter: phi
    [me,~] = calcMeanError(insituData,modelData); % aka bias, greek letter: delta
    [mae,~] = calcMeanAbsoluteError(insituData,modelData);
    mape = calcMeanAbsolutePercentageError(insituData,modelData);

    % Return statistics
    labelMetricNames = {'$r$','slope','bias','RMSE','MAE','MAPE'};
    metrics = [corrCoeff, slope, me, rmse, mae, mape];
    
end % calculatePerformanceMetrics