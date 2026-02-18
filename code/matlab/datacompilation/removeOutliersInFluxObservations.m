function [obsDataNoOutliers,obsErrorNoOutliers] = removeOutliersInFluxObservations(...
    obsData,obsError,nTargetDepths,nLocs)

    % POC flux data are log-normally distributed. When data follow a 
    % log-normal distribution, standard outlier detection methods (like 
    % using the mean and standard deviation directly) are not ideal in the 
    % original scale because of the skewness. Instead, you should take the
    % log10 of your data to transform it into a normal distribution, then 
    % apply standard outlier detection methods.

    obsDataNoOutliers = cellfun(@(x) NaN(size(x)), obsData, 'UniformOutput', false);
    obsErrorNoOutliers = cellfun(@(x) NaN(size(x)), obsError, 'UniformOutput', false);
    nOutliers = 0;
    
    for iTracer = 1:3
        data = obsData{iTracer}; % nTargetDepths x 12 months x nLocs
        error = obsError{iTracer}; 
    
        cleanedData = data;
        cleanedError = error;
        
        for iLoc = 1:nLocs
            for iDh = 1:nTargetDepths

                % Extract time series (e.g., monthly values) for current depth and location
                dataSeries = data(iDh,:,iLoc); 
                errorSeries = error(iDh,:,iLoc);
    
                % Get index of non-NaN values
                isValid = ~isnan(dataSeries);
                validData = dataSeries(isValid); 
    
                % Skip outlier removal if not enough valid data
                if numel(validData) < 2
                    continue;
                end

                % Log-transform the data
                logData = log(validData);
                
                % Compute IQR
                Q1 = prctile(logData, 25);
                Q3 = prctile(logData, 75);
                IQR = Q3 - Q1;
                
                % Define bounds for outliers
                lowerBound = Q1 - 3 * IQR;
                upperBound = Q3 + 3 * IQR;
                
                % Identify outliers
                isOutlier = (logData < lowerBound) | (logData > upperBound);

                % Apply outlier mask back to original vector
                validIdx = find(isValid);  % indices in the original vector that are valid
                outlierIdx = validIdx(isOutlier); % map outlier positions back to original vector
                dataSeries(outlierIdx) = NaN;
                errorSeries(outlierIdx) = NaN;
                
                % Remove outliers from original vectors
                cleanedData(iDh,:,iLoc) = dataSeries;
                cleanedError(iDh,:,iLoc) = errorSeries;
                
                % Print the outlier information: tracer index, depth index, and location index
                if ~isempty(outlierIdx)       
                    fprintf('Outliers removed for Tracer %d, Location %d, Depth %d:\n', iTracer, iLoc, iDh);
                    fprintf('Month indices: %s\n', mat2str(outlierIdx));
                    fprintf('Data before:\n');
                    disp(data(iDh,:,iLoc));
                    fprintf('\nData after:\n');
                    disp(dataSeries);
                    nOutliers = nOutliers + sum(isOutlier);
                end
                
            end % iDh
        end % iLoc
    
        obsDataNoOutliers{iTracer} = cleanedData;
        obsErrorNoOutliers{iTracer} = cleanedError;
    
    end % iTracer

    fprintf('Number of outliers removed: %d\n', nOutliers);

end % removeOutliersInFluxObservations