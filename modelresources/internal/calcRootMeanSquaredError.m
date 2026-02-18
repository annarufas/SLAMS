function [rmse,log_rmse] = calcRootMeanSquaredError(obs,pred)

N = length(obs);

squaredErrors = (pred-obs).^2;
sumOfSquaredErrors = sum(squaredErrors);
rmse = sqrt(sumOfSquaredErrors/N);

squaredLogErrors = (log(pred)-log(obs)).^2; 
sumOfSquaredLogErrors = sum(squaredLogErrors);
log_rmse = sqrt(sumOfSquaredLogErrors/N);

end