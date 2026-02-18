function [mae,log_mae] = calcMeanAbsoluteError(obs,pred)

% Also called mean absolute deviation or difference (MAD) or error (MAE)

N = length(obs);
mae = sum(abs(pred-obs))/N; 
log_mae = sum(abs(log(pred)-log(obs)))/N;

end