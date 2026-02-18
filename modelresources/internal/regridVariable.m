function qData = regridVariable(query_lat,query_lon,query_depth,query_time,...
    source_data,source_lat,source_lon,source_depth,source_time)
    
    % Determine if the depth dimension is present
    if isempty(query_depth) || isempty(source_depth)
        
        % Original grid
        [X, Y, T] = ndgrid(source_lat, source_lon, source_time);
        
        % Query grid
        [qX, qY, qT] = ndgrid(query_lat, query_lon, query_time);

        % Create a gridded interpolant
        F = griddedInterpolant(X, Y, T, source_data, 'linear', 'none');

        % Interpolate to the target grid
        qData = F(qX, qY, qT);
        
    else
        
        % Original grid
        [X, Y, Z, T] = ndgrid(source_lat, source_lon, source_depth, source_time);
        
        % Query grid
        [qX, qY, qZ, qT] = ndgrid(query_lat, query_lon, query_depth, query_time);

        % Create a gridded interpolant
        F = griddedInterpolant(X, Y, Z, T, source_data, 'linear', 'none');

        % Interpolate to the target grid
        qData = F(qX, qY, qZ, qT);
        
    end

end % regridVariable
