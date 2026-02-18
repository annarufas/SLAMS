function y = calculateWeightedAverage(x,nsamples)

% CALCULATEWEIGHTEDAVERAGE Calculates weighted mean.
%
%   INPUT: 
%       x        - observed values
%       nsamples - number of samples that contribute to each observation
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD (2024)
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   This version - Completed 5 April 2024
%
% =========================================================================
%%

totNumSamples = sum(nsamples(:),'omitnan');
accum = 0;
for i = 1:length(x)
    if (~isnan(x(i)))
        accum = accum + x(i).*nsamples(i);
    end
end
y = accum/totNumSamples; 

end