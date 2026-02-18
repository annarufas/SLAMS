function mrpe = calcMeanRelativePercentageError(obs,pred)

N = length(obs);
relativePercentageError = 100.*((pred-obs)./obs);
sumOfRelativePercentageError = sum(relativePercentageError);
mrpe = sumOfRelativePercentageError/N;

end
