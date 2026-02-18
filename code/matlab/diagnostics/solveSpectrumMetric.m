function [metricOutput,metricMCsamples] = solveSpectrumMetric(metricName,...
    spectrumValues_MC,spectrumClasses,nMCsimulations)

% SOLVESPECTRUMMETRIC This function solves for Martin's b, zstar and xi
% using a Monte Carlo eror propagation approach.
% 
%   INPUT:
%       metricName      - String ('b', 'zstar' or 'xi') indicating the metric to solve
%       spectrumValues  - POC flux values (mg C m-2 d-1) or particle number values (# part. m-3 um-1)
%       spectrumClasses - POC flux cast depths (m) or particle number size classes (um)
%       nMCsimulations  - number of Monte Carlo simulations to perform
%
%   OUTPUT:
%       metricOutput    - [median, std upp, std low]
%                         NaN values returned if validation fails
%       metricMCsamples - metric Monte Carlo samples
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
% 
% spectrumValues_MC = dataAvg_MC;
% spectrumClasses = dataClasses;
% nMCsimulations = NUM_MONTE_CARLO_SIMULATIONS;

%% Validate inputs

validMetrics = {'b','zstar','xi'};
if ~ismember(metricName, validMetrics)
    throwError(logID, "Invalid metricName. Choose either 'b', 'zstar' or 'xi'.");
end

%% Define fitting functions

% We need to linearise the power-law and exponential curves by applying
% log10. This way, we can easily propagate the error from a linear
% expression.
fitFuncs = struct(...
    'b',     fittype('-b.*xb + a', 'problem', {'a'}, 'independent', {'xb'}), ...
    'zstar', fittype('-(1/zstar).*xz + a', 'problem', {'a'}, 'independent', {'xz'}), ...
    'xi',    fittype('-xi.*x + a', 'problem', {'a'}, 'independent', {'x'}) ...
);

bounds = struct(...
    'b',     [0.10,    4,  0.09], ...   % [lower, upper, start point]
    'zstar', [  50, 2000,    40], ...
    'xi',    [  -7,    7,  0.10] ...
);

validationParams = struct(...
    'b',     [0.10, 0.10,    4, 10], ...  % [min GOF, lower coeff, upper coeff, max CI factor]
    'zstar', [0.10,   50, 2000, 10], ...
    'xi',    [0.10,   -7,    7, 10] ...
);

%% Initialise output arrays

metric_MC = NaN(nMCsimulations,4); % [value, GOF, upp 95% CI, low 95% CI]
metricOutput = NaN(3,1);           % [median, upp 68% CI, low 68% CI]
metricMCsamples = NaN(nMCsimulations,1);

%% Precompute input data

% The fit function cannot take NaN values, remove them
idxValid = all(~isnan(spectrumValues_MC), 2) & ~isnan(spectrumClasses);
spectrumValues_MC = spectrumValues_MC(idxValid,:);
spectrumClasses = spectrumClasses(idxValid);
nBins = sum(idxValid);

% If not enough data, return NaNs
if nBins < 3
    return;
end

% Precompute log10 values
logSpectrumClasses = log10(spectrumClasses);

% Precompute independent variable based on metric type
x0 = spectrumClasses(1); % Reference size
switch metricName
    case 'b'
        x = logSpectrumClasses - log10(x0);
    case 'zstar'
        x = spectrumClasses - x0;
    otherwise % 'xi'
        x = logSpectrumClasses - log10(x0);
end

% Precompute initial parameters for fitting
a = log10(spectrumValues_MC(1,:));  % a = log10(y0)

%% Start fitting

parfor i = 1:nMCsimulations
    y = log10(spectrumValues_MC(:,i));

    % Fit the data
    [fitResult, gof] = fit(x, y, fitFuncs.(metricName), 'problem', {a(i)}, ...
        'Lower', bounds.(metricName)(1), ...
        'Upper', bounds.(metricName)(2), ...
        'StartPoint', bounds.(metricName)(3));

    % Store fit results
    metric_MC(i,:) = [getfield(fitResult, metricName), gof.adjrsquare, ...
                      max(confint(fitResult)), min(confint(fitResult))];

    % Apply validation 
    metric_MC(i,:) = validateFitOutput(squeeze(metric_MC(i,:)),validationParams.(metricName));
end

%% Store outputs

% Percentiles (16th–84th) give a 68% confidence interval without assuming 
% normality
metricOutput(1) = median(metric_MC(:,1),'omitnan');
metricOutput(2) = prctile(metric_MC(:,1), 84); % upper CI
metricOutput(3) = prctile(metric_MC(:,1), 16); % lower CI
metricMCsamples = metric_MC(:,1);

end % solveSpectrumMetric

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTION
% -------------------------------------------------------------------------

function validated = validateFitOutput(metric,validationParams)

    if isnan(metric(2)) || isinf(metric(2)) || metric(2) < validationParams(1) || ...
       metric(1) < validationParams(2) || metric(1) > validationParams(3) ||...
       metric(3) > validationParams(4) * metric(1)
        validated = NaN(4,1); % return NaN array if validation fails
    else
        validated = metric;
    end
    
end
