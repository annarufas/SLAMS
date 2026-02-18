function generateMetricsLatexOutput(metricScores,locationLabels,varLabels,...
    metricLabels,pathToFile)

    % Reshape misfit data into 2D for LaTeX table
    [nLocs,nVars,nMetrics] = size(metricScores);
    misfitDataReshaped = NaN(nLocs,nVars*nMetrics);
    for iLoc = 1:nLocs
        localMetrics = [];
        for iVar = 1:nVars
            tempMatrix = squeeze(metricScores(iLoc,iVar,:)); 
            localMetrics = [localMetrics;tempMatrix];
        end
        misfitDataReshaped(iLoc,:) = localMetrics;
    end

    % Transpose for LaTeX formatting
    misfitDataTransposed = misfitDataReshaped';

    % Generate the LaTeX table
    dataToOutput = misfitDataTransposed;

    % Open file for writing
    fid = fopen(pathToFile, 'w');
    fprintf(fid, '\\begin{table}[h]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\begin{tabular}{%s}\n', repmat('c', 1, size(dataToOutput, 2)));
    fprintf(fid, '\\toprule\n');

    % Print column headers
    fprintf(fid, '& & %s & %s & %s & %s \\\\\n', locationLabels{:});
    fprintf(fid, '\\midrule\n');

    % Print the data
    iRow = 0;
    for iVar = 1:nVars
        for iMetric = 1:nMetrics

            % Update row
            iRow = iRow + 1;

            % Determine the row data
            rowData = dataToOutput(iRow,:);

            % Format numbers based on value
            formattedRowData = cell(1, numel(rowData));
            for k = 1:numel(rowData)
                if abs(rowData(k)) >= 10
                    formattedRowData{k} = sprintf('%.0f', rowData(k));
                elseif (abs(rowData(k)) < 10 && abs(rowData(k)) >= 1)
                    formattedRowData{k} = sprintf('%.1f', rowData(k));
                else
                    formattedRowData{k} = sprintf('%.2f', rowData(k));
                end
            end

            % If the first statistic in the group
            if mod(iMetric-1, nMetrics) == 0
                % Print a horizontal line before a new group if not the first group
                if iVar > 1
                    fprintf(fid, '\\cmidrule(lr){1-6}\n');
                end
                fprintf(fid, '%s & %s & %s \\\\\n', ...
                    varLabels{iVar}, ...
                    metricLabels{iMetric}, ...
                    sprintf('%s & ', formattedRowData{1:end-1}));
                fprintf(fid, '%s \\\\\n', formattedRowData{end});
            else
                fprintf(fid, '& %s & %s \\\\\n', ...
                    metricLabels{iMetric}, ...
                    sprintf('%s & ', formattedRowData{1:end-1}));
                fprintf(fid, '%s \\\\\n', formattedRowData{end});
            end
        end
    end

    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n');
    fclose(fid);

end % generateMetricsLatexOutput

