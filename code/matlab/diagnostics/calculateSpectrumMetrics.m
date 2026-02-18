function [metricAnnual,metricMonthly] = calculateSpectrumMetrics(metricName,...
    arraySpectrumData,arraySpectrumClasses)

% CALCULATESPECTRUMMETRICS Calculates average Martin's b,
% remineralisation length scale (z*) or the slope of the particle size 
% distribution (xi) and propagates error from input data (POC flux or
% particle size distribution).
%
%   INPUT: 
%       metricName            - string ('b', 'zstar' or 'xi')
%       arraySpectrumData     - spectrum values
%       arraySpectrumClasses  - spectrum classes
%
%   OUTPUT:
%       metricAnnual  - metric, with dimensions: 
%                           for POC flux:         nLocs x 3 (median, uppCI, lowCI) 
%                           for particle numbers: nLocs x 4 x 3 (number of depth horizons // median, uppCI, lowCI) 
%       metricMonthly - metric, with dimensions: 
%                           for POC flux:         nLocs x 12 x 3 (median, uppCI, lowCI) 
%                           for particle numbers: nLocs x 4 x 12 x 3 (number of depth horizons // median, uppCI, lowCI) 
%                              
%   This script uses these external functions:
%       solveSpectrumMetric.m  - custom function
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
        
% To calculate the uncertainty associated to b and z*, it is necessary to
% propagate the error through the fit while considering the uncertainty in 
% the POC flux measurements. Since it's not possible to account for both 
% sources of uncertainty using the 'fittype' and 'worstcase' functions,
% we have opted to perform Monte Carlo sampling across the variable space 
% created by the POC flux error. For each sample, we calculate b and z* 
% and calculate the mean and standard deviation of the ensamble.

%% Definitions

% Parameters for sampling
NUM_MONTE_CARLO_SIMULATIONS = 1e3;
NUM_BOOTSTRAP_SAMPLINGS = 1e5;

% Get number of locations
dims = ndims(arraySpectrumData);
nLocs = size(arraySpectrumData,dims-1);

% Determine output dimensions based on metric type
isXiMetric = strcmp(metricName, 'xi');
depthLevels = 4; % Only relevant for 'xi'

% Define output arrays
if isXiMetric
    metricMonthly = NaN(nLocs,depthLevels,12,3); % 2nd dim: depth horizons, 4th dim: 1=median, 2=68% CI upp, 3=68% CI low
    metricAnnual  = NaN(nLocs,depthLevels,3);
    logResults    = cell(nLocs,depthLevels); 
else
    metricMonthly = NaN(nLocs,12,3); % 3rd dim: 1=median, 2=68% CI upp, 3=68% CI low
    metricAnnual  = NaN(nLocs,3);
    logResults    = cell(nLocs,1);
end 

%% Calculations

% Preallocate local storage for parallel execution
if isXiMetric
    localMetricMonthly = cell(nLocs,depthLevels); 
    localMetricAnnual  = cell(nLocs,depthLevels); 
    localLogResults    = cell(nLocs,1); 
else
    localMetricMonthly = cell(nLocs,1); 
    localMetricAnnual  = cell(nLocs,1); 
    localLogResults    = cell(nLocs,1); 
end

parfor iLoc = 1:nLocs
    fprintf("Worker %d processing iteration %d\n", getCurrentTask().ID, iLoc);
    pause(1);
    
    if isXiMetric
        tempMetricMonthly = NaN(depthLevels,12,3);
        tempMetricAnnual  = NaN(depthLevels,3);
        tempLogResults    = cell(depthLevels,1);
        
        for iDh = 1:depthLevels 
            allMCsamples = NaN(NUM_MONTE_CARLO_SIMULATIONS,12);
            tempMonthly = NaN(12,3);

            for iMonth = 1:12
                % Extract data for current location/month
                dataAvg     = arraySpectrumData(:,iDh,iMonth,iLoc,1); % # part. m-3 um-1
                dataErr     = arraySpectrumData(:,iDh,iMonth,iLoc,2); % # part. m-3 um-1
                dataClasses = arraySpectrumClasses(:,iDh,iMonth,iLoc); % um
                
                if any(~isnan(dataAvg(:)))
                    % Vectorised Monte Carlo sampling
                    dataAvg_MC = dataAvg + dataErr .* randn(size(dataAvg,1), NUM_MONTE_CARLO_SIMULATIONS);
                    dataAvg_MC(dataAvg_MC < 0) = min(dataAvg_MC(dataAvg_MC > 0));  % replace negatives with minimum valid value
                    dataAvg_MC(dataAvg_MC == 0) = NaN; % ensure no zeros before applying log10
    
                    % Calculate metric monthly values
                    [tempMonthly(iMonth,:),allMCsamples(:,iMonth)] =...
                        solveSpectrumMetric(metricName,dataAvg_MC,...
                        dataClasses,NUM_MONTE_CARLO_SIMULATIONS);
                end
            end
            
            % Store temporary results
            tempMetricMonthly(iDh,:,:) = tempMonthly;
            
            % Aggregate monthly results into annual results, use Bootstrapping
            if sum(nnz(tempMonthly(:,1))) >= 2
                bootstrappedSamples = datasample(allMCsamples(:), NUM_BOOTSTRAP_SAMPLINGS, 'Replace', true);

                % Compute confidence intervals
                metricMidval = median(bootstrappedSamples, 'omitnan');
                metricCI = prctile(bootstrappedSamples, [16 84]); % 68% CI
                tempMetricAnnual(iDh,:) = [metricMidval, metricCI(2), metricCI(1)];
                tempLogResults{iDh} = sprintf('\nLocation %d: %s\nEstimated value: %.3f, CI: [%.3f, %.3f]\n', ...
                                     iLoc, metricName, metricMidval, metricCI(1), metricCI(2));
            end

        end % iDh
        
        % Store results in local arrays
        localMetricMonthly{iLoc} = tempMetricMonthly;
        localMetricAnnual{iLoc}  = tempMetricAnnual;
        localLogResults{iLoc}    = tempLogResults;
        
    else
        
        allMCsamples = NaN(NUM_MONTE_CARLO_SIMULATIONS,12);
        tempMonthly = NaN(12,3);  % Temporary array for monthly metrics
        
        for iMonth = 1:12
            % Extract data for current location/month
            dataAvg     = arraySpectrumData(:,iMonth,iLoc,1); % mg C m-2 d-1
            dataErr     = arraySpectrumData(:,iMonth,iLoc,2); % mg C m-2 d-1
            dataClasses = arraySpectrumClasses(:,iMonth,iLoc); % m
            
            if any(~isnan(dataAvg(:)))
                % Vectorised Monte Carlo sampling
                dataAvg_MC = dataAvg + dataErr .* randn(size(dataAvg,1), NUM_MONTE_CARLO_SIMULATIONS);
                dataAvg_MC(dataAvg_MC < 0) = min(dataAvg_MC(dataAvg_MC > 0));  % replace negatives with minimum valid value
                dataAvg_MC(dataAvg_MC == 0) = NaN; % ensure no zeros before applying log10
    
                % Calculate metric monthly values
                [tempMonthly(iMonth,:),allMCsamples(:,iMonth)] = solveSpectrumMetric(...
                    metricName,dataAvg_MC,dataClasses,NUM_MONTE_CARLO_SIMULATIONS);
            end
        end
        
        % Store temporary results
        localMetricMonthly{iLoc} = tempMonthly;
        
        % Aggregate monthly results into annual results, use Bootstrapping
        if sum(nnz(tempMonthly(:,1))) >= 2
            bootstrappedSamples = datasample(allMCsamples(:), NUM_BOOTSTRAP_SAMPLINGS, 'Replace', true);

            % Compute confidence intervals
            metricMidval = median(bootstrappedSamples, 'omitnan');
            metricCI = prctile(bootstrappedSamples, [16 84]); % 68% CI
            localMetricAnnual{iLoc} = [metricMidval, metricCI(2), metricCI(1)];
            localLogResults{iLoc} = sprintf('\nLocation %d: %s\nEstimated value: %.3f, CI: [%.3f, %.3f]\n', ...
                				    iLoc, metricName, metricMidval, metricCI(1), metricCI(2));
        end

    end % metricName
end % iLoc

%% Convert cell arrays back to matrices 
if isXiMetric
    for iLoc = 1:nLocs
        metricMonthly(iLoc,:,:,:) = localMetricMonthly{iLoc};
        metricAnnual(iLoc,:,:) = localMetricAnnual{iLoc};
        logResults(iLoc,:) = localLogResults{iLoc};
    end   
else  
    for iLoc = 1:nLocs
        metricMonthly(iLoc,:,:) = localMetricMonthly{iLoc};
        metricAnnual(iLoc,:) = localMetricAnnual{iLoc};
        logResults{iLoc} = localLogResults{iLoc}; 
    end
end

%% Write logs in one batch
fprintf('%s', strjoin(cellstr(logResults), ''));

end % calculateSpectrumMetrics