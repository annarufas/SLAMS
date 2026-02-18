function [me,log_me] = calcMeanError(obs,pred)

% Also called mean difference (MD), mean error (ME), mean bias (MB) or just
% bias
% Greek symbol: delta

N = length(obs);
me = sum(pred-obs)/N; 
log_me = sum(log(pred)-log(obs))/N;

end