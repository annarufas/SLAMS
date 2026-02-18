function saveFigureInFolder(folderName,figureName)

    % Set default if folderName is not provided or empty
    if nargin < 2 || isempty(folderName)
        targetFolder = fullfile('.','figures');
    else
        targetFolder = fullfile('.','figures',folderName);
    end

    % Create the folder if it doesn't exist
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end

    % Save figure
    exportgraphics(gcf,fullfile(targetFolder,strcat(figureName,'.pdf')),...
        'Resolution',600);
    
    % Optional: PNG format
    exportgraphics(gcf,fullfile(targetFolder,strcat(figureName,'.png')),...
        'Resolution',600);

end