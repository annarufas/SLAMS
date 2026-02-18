function [metricAnnual,metricMonthly] = calculateEfficiencyMetrics(...
    metricName,arrayNumerator,arrayDenominator)

% CALCULATEEFFICIENCYMETRICS Calculates average efficiency metric
% (Teff or PEeff) and propagates error from input data.
%
%   INPUT:
%       metricName       - string ('PEeff' or 'Teff')
%       arrayNumerator   - data array for numerator in quotient, with dimensions: 12 x nLocs x 2 (avg, err)
%       arrayDenominator - data array for denominator in quotient, with same dimensions as above 
%
%   OUTPUT:
%       metricAnnual  - efficiency metric, with dimensions: nLocs x 3 (median, uppCI, lowCI)
%       metricMonthly - efficiency metric, with dimensions: nLocs x 12 x 3 (median, uppCI, lowCI)
%                              
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 14 Feb 2025 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Definitions

% Parameters for sampling
NUM_MONTE_CARLO_SIMULATIONS = 1e3;
NUM_BOOTSTRAP_SAMPLINGS = 1e5;

% Get number of locations
nLocs = size(arrayNumerator,2);

% Define output arrays
metricMonthly = NaN(nLocs,12,3); % 3rd dim: 1=median, 2=68% CI upp, 3=68% CI low
metricAnnual  = NaN(nLocs,3);

%% Start calculations

for iLoc = 1:nLocs
    allMCsamples = NaN(NUM_MONTE_CARLO_SIMULATIONS,12);
    
    % Monthly calculations
    for iMonth = 1:12
        xVals = squeeze(arrayNumerator(iMonth,iLoc,1)); 
        yVals = squeeze(arrayDenominator(iMonth,iLoc,1));
        
        if ~isnan(xVals) && ~isnan(yVals) && (yVals ~= 0)
            exVals = squeeze(arrayNumerator(iMonth,iLoc,2));
            eyVals = squeeze(arrayDenominator(iMonth,iLoc,2));

            % Generate Monte Carlo perturbed samples
            xVals_MC = xVals + exVals .* randn(NUM_MONTE_CARLO_SIMULATIONS, 1);
            yVals_MC = yVals + eyVals .* randn(NUM_MONTE_CARLO_SIMULATIONS, 1);

            % Correct for negative values (set min to smallest positive)
            xVals_MC(xVals_MC < 0) = min(xVals_MC(xVals_MC > 0), [], 'all');
            yVals_MC(yVals_MC < 0) = min(yVals_MC(yVals_MC > 0), [], 'all');

            % Compute the metric for each realisation
            % Identify where both xVals and yVals are not NaN and yVals are not
            % zero (makes the quotient infinite)
            validIndices = ~isnan(xVals_MC) & ~isnan(yVals_MC) & (yVals_MC ~= 0);
            metricMCsamples = NaN(size(xVals_MC));

            % Calculate the quotient for valid entries
            metricMCsamples(validIndices) = xVals_MC(validIndices)./yVals_MC(validIndices);
            allMCsamples(:,iMonth) = metricMCsamples;

            metricMonthly(iLoc,iMonth,1) = median(metricMCsamples(:),'omitnan');
            metricMonthly(iLoc,iMonth,2) = prctile(metricMCsamples(:), 84); % upper CI
            metricMonthly(iLoc,iMonth,3) = prctile(metricMCsamples(:), 16); % lower CI
        end
    end % iMonth
    
    % Aggregate monthly results into annual results, use Bootstrapping
    if sum(nnz(metricMonthly(iLoc,:,1))) >= 2      
        bootstrappedSamples = datasample(allMCsamples(:), NUM_BOOTSTRAP_SAMPLINGS, 'Replace', true);

        % Compute confidence intervals
        metricMidval = median(bootstrappedSamples, 'omitnan');
        metricCI = prctile(bootstrappedSamples, [16 84]); % 68% CI
        metricAnnual(iLoc,:) = [metricMidval, metricCI(2), metricCI(1)];

        fprintf('\nLocation %d: %s\nEstimated value: %.3f, CI: [%.3f, %.3f]\n', ...
            iLoc, metricName, metricMidval, metricCI(1), metricCI(2));
    end

end % iLoc
    
end % calculateEfficiencyMetrics