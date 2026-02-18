function dataConverted = applyDensityConversion(data,rho,data_lat,data_lon,...
    data_depth,rho_lat,rho_lon,rho_depth,t)

    % Input validation
    if ndims(data) ~= 4 && (size(data,1) > size(data,2))
        error('Input data array must be a 4D matrix of dimensions lat x lon x depth x time.');
    end
    if ndims(rho) ~= 4 && (size(rho,1) > size(rho,2))
        error('Input rho array must be a 4D matrix of dimensions lat x lon x depth x time.');
    end

    % Regrid rho to same grid used by data
    if ~isequal(size(data), size(rho))
        rho = regridVariable(data_lat,data_lon,data_depth,t,...
            rho,rho_lat,rho_lon,rho_depth,t); % kg m-3
    end

    % Convert
    dataConverted = data.*rho.*1e-3; % umol kg-1 --> umol L-1 (=mmol m-3)

end % applyDensityConversion