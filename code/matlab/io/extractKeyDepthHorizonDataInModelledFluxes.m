function [fluxMonthlyDh,fluxMonthlyDhDepths,fluxAnnualDh,fluxAnnualDhDepths] =...
    extractKeyDepthHorizonDataInModelledFluxes(...
        config,fluxMonthlyProfile,keyDepthLevelsMonthlyFluxData)

    dims = size(fluxMonthlyProfile);
    nDims = ndims(fluxMonthlyProfile);

    nSampledDepths = dims(1);
    nMonths = dims(2);
    nLocs = dims(3);

    % Determine whether the 4th dimension is 'stat' (mean/std) or tracer
    if nDims == 4
        % Case: (depth, month, loc, tracer)
        hasStats = false;
        nTracers = dims(4);
    elseif nDims == 5
        % Case: (depth, month, loc, stat, tracer)
        hasStats = true;
        nStats = dims(4);  % Expecting nStats = 2
        nTracers = dims(5);
    else
        error('Unexpected dimensions in seqFluxMonthlyProfile');
    end

    nTargetDepths = size(keyDepthLevelsMonthlyFluxData,1);

    fluxMonthlyDh = NaN(nTargetDepths,nMonths,nLocs,2,nTracers); % 1st dim: 1=zeu, 2=zmidmeso, 3=zmeso, 4=zbathy, 4th dim: 1=avg, 2=err
    fluxMonthlyDhDepths = NaN(nTargetDepths,nMonths,nLocs,nTracers);
    fluxAnnualDh = NaN(nTargetDepths,nLocs,2,nTracers); 
    fluxAnnualDhDepths = NaN(nTargetDepths,nLocs,nTracers);

    validDepthMask = false(nSampledDepths,nMonths,nLocs,nTargetDepths,nTracers);
    for iTracer = 1:nTracers
        for iDh = 1:nTargetDepths
            for iMonth = 1:nMonths
                for iLoc = 1:nLocs

                    if (ndims(keyDepthLevelsMonthlyFluxData) == 5) % contains tracer-specific depths
                        depthRangeTop = keyDepthLevelsMonthlyFluxData(iDh,iMonth,iLoc,1,iTracer);    % upper boundary depth
                        depthRangeBottom = keyDepthLevelsMonthlyFluxData(iDh,iMonth,iLoc,2,iTracer); % lower boundary depth
                    elseif (ndims(keyDepthLevelsMonthlyFluxData) == 4) % does not contain tracer-specific depths
                        depthRangeTop = keyDepthLevelsMonthlyFluxData(iDh,iMonth,iLoc,1);    % upper boundary depth
                        depthRangeBottom = keyDepthLevelsMonthlyFluxData(iDh,iMonth,iLoc,2); % lower boundary depth
                    end

                    % Only proceed if both are numbers ratehr than NaN
                    if ~isnan(depthRangeTop) && ~isnan(depthRangeBottom)
                        % Find the closest depth greater than or equal to depthRangeTop
                        [~,idxTop] = min(abs(config.availSedTrapDeployDepths - depthRangeTop));
        
                        % Find the closest depth less than or equal to depthRangeBottom
                        [~,idxBottom] = min(abs(config.availSedTrapDeployDepths - depthRangeBottom));
                        if idxBottom > nSampledDepths
                            idxBottom = nSampledDepths;
                        end
        
                        % Extract valid modelled depth range
                        validDepthMask((idxTop:idxBottom),iMonth,iLoc,iDh,iTracer) = true;
        
                        % Calculate the average monthly unique modelled depth
                        fluxMonthlyDhDepths(iDh,iMonth,iLoc,iTracer) =...
                            (config.availSedTrapDeployDepths(idxTop)+config.availSedTrapDeployDepths(idxBottom))/2;
                    end
                end % iLoc
            end % iMonth
    
            % Calculate the average annual unique modelled depth
            fluxAnnualDhDepths(iDh,iLoc,iTracer) = mean(fluxMonthlyDhDepths(iDh,:,iLoc,iTracer),2,'omitnan');
        end % iDh
    end % iTracer

    % Monthly calculations, with error propagation
    for iTracer = 1:nTracers
        for iDh = 1:nTargetDepths
            % Mask for the current depth layer (size nSampledDepths x nMonths x nLocs)
            currentMask = squeeze(validDepthMask(:,:,:,iDh,iTracer));

            if hasStats
                fluxValues = squeeze(fluxMonthlyProfile(:,:,:,1,iTracer)); % mean
                fluxErrors = squeeze(fluxMonthlyProfile(:,:,:,2,iTracer)); % fractional std
                fluxErrors = fluxErrors .* fluxValues; % compute actual error values (currently those are fractional values)
            else
                fluxValues = squeeze(fluxMonthlyProfile(:,:,:,iTracer)); % only mean, no std
                fluxErrors = NaN(size(fluxValues)); % no errors available
            end

            % Select only valid positions based on the precomputed mask
            fluxValues(~currentMask) = NaN;
            fluxErrors(~currentMask) = NaN;

            % Compute mean flux excluding NaNs
            fluxMonthlyDh(iDh,:,:,1,iTracer) = squeeze(mean(fluxValues, 1, 'omitnan'));

            % Store propagated error if available
            if hasStats
                for iMonth = 1:nMonths    
                    validFlux = squeeze(fluxValues(:,iMonth,:)); % nSampledDepths x nLocs
                    validError = squeeze(fluxErrors(:,iMonth,:));
                    if any(~isnan(validError(:)))
                        fluxMonthlyDh(iDh,iMonth,:,2,iTracer) = ...
                            calculateUncertaintyInAverage(validFlux,validError);
                    end
                end
            end
        end % iDh
    end % iTracer

    % Annual average, with error propagation
    for iTracer = 1:nTracers
        for iDh = 1:nTargetDepths
            fluxValues = squeeze(fluxMonthlyDh(iDh,:,:,1,iTracer));
            fluxAnnualDh(iDh,:,1,iTracer) = squeeze(mean(fluxValues, 1, 'omitnan')); % exclude NaNs across months

            if hasStats
                fluxErrors = squeeze(fluxMonthlyDh(iDh,:,:,2,iTracer));
                if any(~isnan(fluxErrors(:)))
                    fluxAnnualDh(iDh,:,2,iTracer) = ...
                        calculateUncertaintyInAverage(fluxValues, fluxErrors);
                end
            end
        end
    end

end % extractKeyDepthHorizonDataInModelledFluxes