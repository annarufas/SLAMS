function [pnumMonthlyDh,pnumMonthlyDhDepths,pnumAnnualDh,pnumAnnualDhDepths] =...
    extractKeyDepthHorizonDataInModelledParticleNumbers(config,...
        pnumMonthlyProfile,keyDepthLevelsMonthlyData)

    dims = size(pnumMonthlyProfile);
    nDims = ndims(pnumMonthlyProfile);

    nSizeClasses = dims(1);
    nSampledDepths = dims(2);
    nMonths = dims(3);
    nLocs = dims(4);

    if nDims == 4
        hasStats = false;
    elseif nDims == 5
        hasStats = true;
    else
        error('Unexpected dimensions in seqPnumMonthlyProfile');
    end

    nTargetDepths = size(keyDepthLevelsMonthlyData,1);

    pnumMonthlyDh = NaN(nSizeClasses,nTargetDepths,nMonths,nLocs,2); % 2nd dim: 1=zeu, 2=zmidmeso, 3=zmeso, 4=zbathy, 5th dim: 1=avg, 2=err
    pnumMonthlyDhDepths = NaN(nTargetDepths,nMonths,nLocs); 
    pnumAnnualDh = NaN(nSizeClasses,nTargetDepths,nLocs,2); 
    pnumAnnualDhDepths = NaN(nTargetDepths,nLocs); 
    
    validDepthMask = false(nSampledDepths,nMonths,nLocs,nTargetDepths);
    for iDh = 1:nTargetDepths
        for iMonth = 1:nMonths
            for iLoc = 1:nLocs
                depthRangeTop = keyDepthLevelsMonthlyData(iDh,iMonth,iLoc,1);    % upper boundary depth
                depthRangeBottom = keyDepthLevelsMonthlyData(iDh,iMonth,iLoc,2); % lower boundary depth

                % Find the closest depth greater than or equal to depthRangeTop
                [~,idxTop] = min(abs(config.availImagSysDeployDepths - depthRangeTop));

                % Find the closest depth less than or equal to depthRangeBottom
                [~,idxBottom] = min(abs(config.availImagSysDeployDepths - depthRangeBottom));
                if idxBottom > nSampledDepths
                    idxBottom = nSampledDepths;
                end
    
                % Extract valid modelled depth range
                validDepthMask((idxTop:idxBottom),iMonth,iLoc,iDh) = true;

                % Calculate the average monthly unique modelled depth
                pnumMonthlyDhDepths(iDh,iMonth,iLoc) =...
                    (config.availImagSysDeployDepths(idxTop)+config.availImagSysDeployDepths(idxBottom))/2;
            end % iLoc
        end % iMonth

        % Calculate the average annual unique modelled depth
        pnumAnnualDhDepths(iDh,iLoc) = mean(pnumMonthlyDhDepths(iDh,:,iLoc),2,'omitnan');
    end
    
    % Monthly calculations, with error propagation
    for iSc = 1:nSizeClasses
        for iDh = 1:nTargetDepths
            % Mask for the current depth layer (size nSampledDepths x nMonths x nLocs)
            currentMask = squeeze(validDepthMask(:,:,:,iDh));

            if hasStats
                pnumValues = squeeze(pnumMonthlyProfile(iSc,:,:,:,1));
                pnumErrors = squeeze(pnumMonthlyProfile(iSc,:,:,:,2));
                pnumErrors = pnumErrors .* pnumValues; % compute actual error values (currently those are fractional values)
            else
                pnumValues = squeeze(pnumMonthlyProfile(iSc,:,:,:));
                pnumErrors = NaN(size(pnumValues));
            end

            % Select only valid positions based on the precomputed mask
            pnumValues(~currentMask) = NaN;
            pnumErrors(~currentMask) = NaN;

            % Compute mean particle number excluding NaNs
            pnumMonthlyDh(iSc,iDh,:,:,1) = squeeze(mean(pnumValues, 1, 'omitnan'));

            % Error propagation
            if hasStats
                for iMonth = 1:nMonths     
                    validPnum = squeeze(pnumValues(:,iMonth,:)); % nSampledDepths x nLoc
                    validError = squeeze(pnumErrors(:,iMonth,:));
                    if any(~isnan(validError(:)))
                        pnumMonthlyDh(iSc,iDh,iMonth,:,2) = ...
                            calculateUncertaintyInAverage(validPnum,validError);
                    end
                end
            end
        end % iDh
    end % iSc

    % Annual average, with error propagation
    for iSc = 1:nSizeClasses
        for iDh = 1:nTargetDepths
            pnumValues = squeeze(pnumMonthlyDh(iSc,iDh,:,:,1));
            pnumAnnualDh(iSc,iDh,:,1) = squeeze(mean(pnumValues, 1, 'omitnan')); % compute mean flux excluding NaNs across months
            
            if hasStats
                pnumErrors = squeeze(pnumMonthlyDh(iSc,iDh,:,:,2));
                if any(~isnan(pnumErrors(:)))
                    pnumAnnualDh(iSc,iDh,:,2) = ...
                        calculateUncertaintyInAverage(pnumValues,pnumErrors);
                end
            end
        end
    end

end % extractKeyDepthHorizonDataInModelledParticleNumbers