function [outputData,caxisMin,caxisMax,ticks,tickLabels] =...
    assessLogarithmicTransformation(originalData)

% originalData = data;

    % This function:
    % - handling of logarithmic and linear scaling separately
    % - dynamically adjusts step sizes to ensure a maximum of 8 ticks
    % - handles small ranges by refining limits and step sizes
    % - robust percentile handling to prevent extreme outliers from skewing the limits

    useLogScale = false;
    
    originalData(originalData==0) = NaN; % make sure there are no zeros(log10 cannot cope with those)

    % Compute log range and skewness
    logRange = log10(max(originalData(:))) - log10(min(originalData(:)));
    dataSkewness = skewness(originalData(:));
    % figure()    
    % histogram(originalData, 20);     % Plot histogram with 20 bins
    % xlabel('Data values');
    % ylabel('Frequency');
    % title('Histogram of Data')

    % Use log scale if data spans at least 4 orders of magnitude OR is highly skewed
    % if logRange > 3 || dataSkewness > 1
    if logRange > 6 || dataSkewness > 1
        useLogScale = true;
    end

    % **LOGARITHMIC SCALE HANDLING**
    if useLogScale
        outputData = log10(originalData);
        logData = log10(originalData(:));

        % Compute 5th and 95th percentiles to exclude extreme outliers
        logMin = prctile(logData, 5);
        logMax = prctile(logData, 95);
          
        % Round to nearest power of 10
        caxisMin = floor(logMin);
        caxisMax = floor(logMax);

        % Adjust max limit slightly if needed
        if (max(logData) - caxisMax) < 0.5
            caxisMax = caxisMax + 1;
        end

        % Determine step size to ensure max 8 ticks
        nTicks = caxisMax - caxisMin + 1;
        if nTicks > 8
            stepSize = ceil((caxisMax - caxisMin) / 8);
        else
            stepSize = 1;
        end

        % Generate tick positions
        ticks = caxisMin:stepSize:caxisMax;

        % Create labels
        tickLabels = arrayfun(@(x) sprintf('10^{%d}',x), ticks, 'UniformOutput', false);
        
        %%
%         % Adjust min and max to nearest power of 10
%         caxisMin = 10^floor(logMin);
%         caxisMax = 10^ceil(logMax);
% 
%         % Ensure max > min
%         if caxisMax <= caxisMin
%             caxisMax = caxisMin * 10;
%         end
% 
%         % Determine step size dynamically (max 8 ticks)
%         nTicks = 8;
%         tickStep = (log10(caxisMax) - log10(caxisMin)) / (nTicks - 1);
%         ticks = log10(caxisMin):tickStep:log10(caxisMax);
%         ticks = 10.^ticks; % Convert back to normal scale
% 
%         % Create labels in 10^x format
%         tickLabels = arrayfun(@(x) sprintf('10^{%d}', round(log10(x))), ticks, 'UniformOutput', false);
        %%
    else
        % **LINEAR SCALE HANDLING**
        outputData = originalData;

        % Compute 5th and 95th percentiles for robust min/max
        linearMin = prctile(originalData(:), 5);
        linearMax = prctile(originalData(:), 95);
        
        % Ensure min and max are not identical
        if linearMax == linearMin
            linearMax = linearMax * 1.1;
        end
        
        % **Handle very small or very large numbers**
        range = linearMax - linearMin;
        if range < 1e-6
            
            % Define caxisMin and caxisMax properly in linear space
            caxisMin = floor(linearMin * 10^ceil(-log10(linearMax))) / 10^ceil(-log10(linearMax));
            caxisMax = ceil(linearMax * 10^ceil(-log10(linearMax))) / 10^ceil(-log10(linearMax));

            % Ensure at least 3 ticks
            nTicks = max(3, min(8, ceil(log10(caxisMax) - log10(caxisMin)) + 1));

            % Generate tick positions with appropriate spacing
            ticks = linspace(caxisMin, caxisMax, nTicks);

            % Format tick labels using scientific notation
            tickLabels = arrayfun(@(x) sprintf('%.1e', x), ticks, 'UniformOutput', false);
    
        else
    
            % Use standard rounding for normal ranges
            roundFactor = 10^floor(log10(range));
            caxisMin = floor(linearMin / roundFactor) * roundFactor;
            caxisMax = ceil(linearMax / roundFactor) * roundFactor;

            % Ensure max > min
            if caxisMax <= caxisMin
                caxisMax = caxisMin + roundFactor;
            end

            % Determine step size dynamically (max 8 ticks)
            stepSizes = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
            stepSize = stepSizes(find((caxisMax - caxisMin) ./ stepSizes <= 8, 1, 'first'));

            % Generate tick positions with at least 3 ticks but no more
            % than 8
            ticks = caxisMin:stepSize:caxisMax;
            if numel(ticks) < 3
                ticks = linspace(caxisMin, caxisMax, 3);
            elseif numel(ticks) > 8
                ticks = linspace(caxisMin, caxisMax, 8);
            end
            
            % Format tick labels
            % Ensure caxisMin is not zero (log10(0) is undefined)
            if caxisMin == 0
                caxisMin = min(ticks(ticks > 0)); % Set to the smallest positive tick
            end
            if abs(log10(caxisMin)) >= 4 || abs(log10(caxisMax)) >= 4
                tickLabels = arrayfun(@(x) sprintf('%.1e', x), ticks, 'UniformOutput', false);
            else
                tickLabels = arrayfun(@(x) sprintf('%g', x), ticks, 'UniformOutput', false);
            end

        end

    end

end % assessLogarithmicTransformation