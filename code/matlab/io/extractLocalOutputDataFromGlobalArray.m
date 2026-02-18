function subsetOutput = extractLocalOutputDataFromGlobalArray(config,nLocs,...
    listLocalLons,listLocalLats,output)

    % If we have run the model globally, we can extract local data of interest

    % Output index for each target point
    selectedIdxs = NaN(nLocs,1);

    for iLoc = 1:nLocs

        % Find the idx to the coordinates in the global grid
        [~,iLon] = min(abs(config.lons(:) - listLocalLons(iLoc)));
        [~,iLat] = min(abs(config.lats(:) - listLocalLats(iLoc)));

        lonToSearch = config.lons(iLon);
        latToSearch = config.lats(iLat);

        isMatch = (config.gridLons == lonToSearch) & (config.gridLats == latToSearch);
        idxLoc = find(isMatch);
        
        selectedIdxs(iLoc) = idxLoc;
        
    end    

    % Get field names
    fields = fieldnames(output);
    
    % Initialise new structure to hold extracted data
    subsetOutput = struct();
    
    % Loop through each field
    for f = 1:numel(fields)
        fieldName = fields{f};
        data = output.(fieldName);
    
        % Get number of dimensions
        dims = ndims(data);
    
        % Build indexing: ':' for all other dimensions, and idxClosest for the last one
        idx = repmat({':'}, 1, dims);  % e.g., {':', ':', ':'}
        idx{dims} = selectedIdxs;      % replace last index with your location indices
    
        % Extract subset
        subsetOutput.(fieldName) = data(idx{:});
    end

end % extractLocalOutputDataFromGlobalArray