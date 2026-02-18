function e = calculateUncertaintyInWeightedAverage(xVals,eVals,nVals)

N = length(xVals);

nTot = sum(nVals(:),'omitnan'); 
varianceWeighted = NaN(N,1);
for iVal = 1:N 
    if (~isnan(xVals(iVal)))
        normWeight = nVals(iVal)/nTot; % normalised weight
        varianceWeighted(iVal) = normWeight^2 * eVals(iVal)^2;
    end
end

e = sqrt(sum(varianceWeighted(:),'omitnan'));

end
