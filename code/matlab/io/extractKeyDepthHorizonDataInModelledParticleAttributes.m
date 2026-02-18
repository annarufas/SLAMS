function pattAnnualDh = extractKeyDepthHorizonDataInModelledParticleAttributes(...
    config,pattAnnualProfile,keyDepthLevelsAnnualData)

    [nSizeCategories,nSampledDepths,nAttributes,nLocs] = size(pattAnnualProfile);
    nTargetDepths = size(keyDepthLevelsAnnualData,1);
    pattAnnualDh = NaN(nSizeCategories,nTargetDepths,nAttributes,nLocs); % 2nd dim: 1=zeu, 2=zmidmeso, 3=zmeso, 4=zbathy
    
    validDepthMask = false(nSampledDepths,nLocs,nTargetDepths);
    for iDh = 1:nTargetDepths
        for iLoc = 1:nLocs
            depthRangeTop = keyDepthLevelsAnnualData(iDh,iLoc,1);    % upper boundary depth
            depthRangeBottom = keyDepthLevelsAnnualData(iDh,iLoc,2); % lower boundary depth

            % Find the closest depth greater than or equal to depthRangeTop
            [~,idxTop] = min(abs(config.availImagSysDeployDepths - depthRangeTop));

            % Find the closest depth less than or equal to depthRangeBottom
            [~,idxBottom] = min(abs(config.availImagSysDeployDepths - depthRangeBottom));
            if idxBottom > nSampledDepths
                idxBottom = nSampledDepths;
            end
    
            % Extract valid modelled depth range
            validDepthMask((idxTop:idxBottom),iLoc,iDh) = true;
        end
    end
    
    % Extract
    for iSc = 1:nSizeCategories
        for iDh = 1:nTargetDepths
            for iAtt = 1:nAttributes
                
                % Mask for the current depth layer (size nSampledDepths x nLocs)
                currentMask = squeeze(validDepthMask(:,:,iDh));
                pattValues = squeeze(pattAnnualProfile(iSc,:,iAtt,:));

                % Select only valid positions based on the precomputed mask
                pattValues(~currentMask) = NaN;

                % Compute mean particle number excluding NaNs
                pattAnnualDh(iSc,iDh,iAtt,:) = squeeze(mean(pattValues, 1, 'omitnan'));

            end % iAtt
        end % iDh
    end % iSc
    
end % extractKeyDepthHorizonDataInModelledParticleAttributes
