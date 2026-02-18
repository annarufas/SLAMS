function [fluxRelUncertainty, pNumRelUncertainty] = calculateBootstrapUncertainty(seqOutput)
% Outputs:
%   - totalRelUncertainty : relative uncertainty (fractional, 0–1)
%                           NOT percent (0–100)

    nBootstraps = 1000; 
    nBootstrapSimulations = 50;

    fluxRelUncertainty = NaN(3, nBootstrapSimulations);
    pNumRelUncertainty = NaN(nBootstrapSimulations, 1);

    for i = 1:nBootstrapSimulations
        [~, fluxRelUncertainty(:, i)] = bootstrapFluxUncertainty(seqOutput, nBootstraps);
        [~, pNumRelUncertainty(i)] = bootstrapParticleNumberUncertainty(seqOutput, nBootstraps);
    end

    fluxRelUncertainty = mean(fluxRelUncertainty, 2, 'omitnan');
    pNumRelUncertainty = mean(pNumRelUncertainty, 'omitnan');

end

function [bootResults,totalRelUncertainty] = bootstrapFluxUncertainty(...
   seqOutput,nBootstraps)

    % Bootstrapping on model outputs to estimate uncertainty in particle 
    % flux data without rerunning the full model.
    % Inputs:
    %   - seqOutput: structure with monthly flux data
    %   - nBootstraps: number of bootstrap iterations
    % Outputs:
    %   - bootResults: Struct containing bootstrapped means, confidence intervals, and relative uncertainties
    %   - totalRelUncertainty: Total combined uncertainty across dimensions
    %     (relative percentage, 0-1)

    % Define tracers
    tracers = struct(...
        'POC', seqOutput.monthlyOrgCarbonFlux + seqOutput.monthlyTepFlux, ...
        'PIC', seqOutput.monthlyCalciteFlux, ...
        'BSi', seqOutput.monthlyOpalFlux);
    tracerNames = fieldnames(tracers);
    nTracers = numel(tracerNames);

    % Initialise result struct
    bootResults = struct();
    dims = {'Spatial', 'Temporal', 'Depth'};

    for iTracer = 1:nTracers
        fluxData = tracers.(tracerNames{iTracer});
        fluxData(fluxData == 0) = NaN; % filter out zeros sicne they represent no data values
        results = struct();

        % Compute bootstrap statistics for spatial, temporal, and depth
        results.Spatial  = bootstrapDimension(mean(reshape(fluxData, [], size(fluxData, 3)), 1, 'omitnan'), nBootstraps);
        results.Temporal = bootstrapDimension(squeeze(mean(fluxData, [1, 3], 'omitnan')), nBootstraps);
        results.Depth    = bootstrapDimension(squeeze(mean(fluxData, [2, 3], 'omitnan')), nBootstraps);
 
        % Compute total relative uncertainty
        ruValues = cellfun(@(dim) results.(dim).RelUncertainty, dims);
        results.TotalRelUncertainty = sqrt(sum(ruValues.^2));

        % Store in results struct
        bootResults.(tracerNames{iTracer}) = results;
    end
    
    % Extract total uncertainties for all tracers
    totalRelUncertainty = cellfun(@(tr) bootResults.(tr).TotalRelUncertainty, tracerNames);

    function result = bootstrapDimension(data, nBootstraps)
        % Bootstrapping function for a given data dimension
        data = data(~isnan(data));  % remove NaNs
        bootMeans = arrayfun(@(~) mean(datasample(data, length(data), 'Replace', true)), 1:nBootstraps);
        ci = prctile(bootMeans, [16 84]);
        result.Mean = mean(bootMeans);
        result.CI = ci;
        result.RelUncertainty = (ci(2) - ci(1)) / (2 * result.Mean);
    end

end % bootstrapFluxUncertainty

function [bootResults,totalRelUncertainty] = bootstrapParticleNumberUncertainty(...
    seqOutput,nBootstraps)

    % Function to estimate uncertainty in particle numbers using bootstrapping
    % Inputs:
    %   - particleNumbers: 16 x 4 x 12 x 4448 array (size class x depth x month x location)
    %   - nBootstraps: number of bootstrap iterations
    % Outputs:
    %   - bootResults: struct containing bootstrapped means & confidence intervals
    %   - totalRelUncertainty: combined uncertainty across all factors
    %     (percentage, 0-1)

    % Array
    particleNumbers = seqOutput.monthlyParticleNumInSc;

    % Initialise storage for results
    bootResults = struct();

    % Spatial bootstrapping
    spatialMeans = zeros(nBootstraps,1);
    locationMeans = squeeze(mean(particleNumbers, [1,2,3], 'omitnan')); % average over size classes, depths & months
    locationMeans = locationMeans(~isnan(locationMeans));  % remove NaNs
    for i = 1:nBootstraps
        resampledLocs = datasample(locationMeans, size(locationMeans,2), 'Replace', true);
        spatialMeans(i) = mean(resampledLocs, 'omitnan');
    end
    spatialCI = prctile(spatialMeans, [16 84]);
    bootResults.Spatial = struct('Mean', mean(spatialMeans), ...
                                 'CI', spatialCI, ...
                                 'RelativeUncertainty', (spatialCI(2) - spatialCI(1)) / (2 * mean(spatialMeans)));

    % Temporal bootstrapping
    temporalMeans = zeros(nBootstraps,1);
    monthlyMeans = squeeze(mean(particleNumbers, [1,2,4], 'omitnan')); % average over sizes, depths and locations
    monthlyMeans = monthlyMeans(~isnan(monthlyMeans));  % remove NaNs
    for i = 1:nBootstraps
        resampledMonths = datasample(monthlyMeans, length(monthlyMeans), 'Replace', true);
        temporalMeans(i) = mean(resampledMonths, 'omitnan');
    end
    temporalCI = prctile(temporalMeans, [16 84]);
    bootResults.Temporal = struct('Mean', mean(temporalMeans), ...
                                  'CI', temporalCI, ...
                                  'RelativeUncertainty', (temporalCI(2) - temporalCI(1)) / (2 * mean(temporalMeans)));

    % Depth bootstrapping
    depthMeans = zeros(nBootstraps,1);
    depthLevels = squeeze(mean(particleNumbers, [1,3,4], 'omitnan')); % average over size classes, months and locations
    depthLevels = depthLevels(~isnan(depthLevels));  % remove NaNs
    for i = 1:nBootstraps
        resampledDepths = datasample(depthLevels, length(depthLevels), 'Replace', true);
        depthMeans(i) = mean(resampledDepths, 'omitnan');
    end
    depthCI = prctile(depthMeans, [16 84]);
    bootResults.Depth = struct('Mean', mean(depthMeans), ...
                               'CI', depthCI, ...
                               'RelativeUncertainty', (depthCI(2) - depthCI(1)) / (2 * mean(depthMeans)));

    % Size class bootstrapping
    sizeClassMeans = zeros(nBootstraps,1);
    sizeClassLevels = squeeze(mean(particleNumbers, [2,3,4], 'omitnan')); % average over depths, months & locations
    sizeClassLevels = sizeClassLevels(~isnan(sizeClassLevels));  % remove NaNs
    for i = 1:nBootstraps
        resampledSizes = datasample(sizeClassLevels, length(sizeClassLevels), 'Replace', true);
        sizeClassMeans(i) = mean(resampledSizes, 'omitnan');
    end
    sizeClassCI = prctile(sizeClassMeans, [16 84]);
    bootResults.SizeClass = struct('Mean', mean(sizeClassMeans), ...
                                   'CI', sizeClassCI, ...
                                   'RelativeUncertainty', (sizeClassCI(2) - sizeClassCI(1)) / (2 * mean(sizeClassMeans)));

    % Combine uncertainties
    totalRelUncertainty = sqrt(...
        bootResults.Spatial.RelativeUncertainty^2 + ...
        bootResults.Temporal.RelativeUncertainty^2 + ...
        bootResults.Depth.RelativeUncertainty^2 + ...
        bootResults.SizeClass.RelativeUncertainty^2);

end % bootstrapParticleNumberUncertainty
