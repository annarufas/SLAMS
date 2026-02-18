function mape = calcMeanAbsolutePercentageError(obs,pred)

% Also called mean absolute percentage difference (MAPD) or mean absolute
% percentage error (MAPE)

N = length(obs);
absolutePercentageError = 100.*(abs(pred-obs)./obs);
sumOfAbsolutePercentageError = sum(absolutePercentageError);
mape = sumOfAbsolutePercentageError/N;

end
