function integratedData = integrateOverDepth(data,config)

    depthLimit = 200; % m, depth limit for integration
    dataSize = size(data);
    hasMonthlyDimension = (length(dataSize) == 4);  % check if the data is 4D (i.e., monthly array instead of annual)
    integratedData = NaN(config.nLats,config.nLons,hasMonthlyDimension*12+(1-hasMonthlyDimension)); % mg C m-2
    
    for iLon = 1:config.nLons
        for iLat = 1:config.nLats
            localDepths = squeeze(config.geoDepths(iLat,iLon,:));
            
            if hasMonthlyDimension
                % Handle 4D data (lat x lon x depth x month)
                for iMonth = 1:12
                    localConc = squeeze(data(iLat,iLon,:,iMonth)); % mg C m-3
                    validDepths = localDepths <= depthLimit;
                    
                    if any(validDepths) && any(localConc(validDepths) > 0)
                        % Restrict to valid depths and concentrations
                        depthsToIntegrate = localDepths(validDepths);
                        concToIntegrate = localConc(validDepths);

                        % Perform trapezoidal integration
                        integratedData(iLat,iLon,iMonth) = trapz(depthsToIntegrate,concToIntegrate); % mg C m-2
                        if isnan(integratedData(iLat,iLon,iMonth))
                            throwError('ERROR: integration result is NaN for month %d at lat %d, lon %d.\n', iMonth, iLat, iLon);
                        end
                    end
                end
                
            else
                % Handle 3D data (lat x lon x depth)
                localConc = squeeze(data(iLat,iLon,:)); % mg C m-3
                validDepths = localDepths <= depthLimit;

                if any(validDepths) && any(localConc(validDepths) > 0)
                    % Restrict to valid depths and concentrations
                    depthsToIntegrate = localDepths(validDepths);
                    concToIntegrate = localConc(validDepths);

                    % Perform trapezoidal integration
                    integratedData(iLat,iLon) = trapz(depthsToIntegrate,concToIntegrate); % mg C m-2
                    if isnan(integratedData(iLat,iLon))
                        throwError('ERROR: integration result is NaN for lat %d, lon %d.\n.', iLat, iLon);
                    end
                end
            end
        end % iLat
    end % iLon

end % integrateOverDepth