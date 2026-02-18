function [smsAnnualDh,auxAnnualDh] = extractKeyDepthHorizonDataInModelledSmsAndAuxTerms(...
    config,smsAnnualProfile,auxAnnualProfile,keyDepthLevelsAnnualData)

    [nSmsTerms,nSampledDepths,nLocs] = size(smsAnnualProfile);
    [nAuxTerms,~,~] = size(auxAnnualProfile);
    nTargetDepths = size(keyDepthLevelsAnnualData,1);
    
    smsAnnualDh = NaN(nSmsTerms,nTargetDepths,nLocs); % 2nd dim: 1=zeu, 2=zmidmeso, 3=zmeso, 4=zbathy
    auxAnnualDh = NaN(nAuxTerms,nTargetDepths,nLocs); % 2nd dim: 1=zeu, 2=zmidmeso, 3=zmeso, 4=zbathy

    validDepthMask = false(nSampledDepths,nLocs,nTargetDepths);
    for iDh = 1:nTargetDepths
        for iLoc = 1:nLocs
            depthRangeTop = keyDepthLevelsAnnualData(iDh,iLoc,1);    % upper boundary depth
            depthRangeBottom = keyDepthLevelsAnnualData(iDh,iLoc,2); % lower boundary depth

            % Find the closest depth greater than or equal to depthRangeTop
            [~,idxTop] = min(abs(config.gridDepths(:,iLoc) - depthRangeTop));

            % Find the closest depth less than or equal to depthRangeBottom
            [~,idxBottom] = min(abs(config.gridDepths(:,iLoc) - depthRangeBottom));

            % Extract valid modelled depth range
            validDepthMask((idxTop:idxBottom),iLoc,iDh) = true;
        end
    end
    
    % Process SMS terms
    for iTerm = 1:nSmsTerms
        for iDh = 1:nTargetDepths
            % Apply precomputed mask
            currentMask = squeeze(validDepthMask(:,:,iDh));
            smsValues = squeeze(smsAnnualProfile(iTerm,:,:));
            smsValues(~currentMask) = NaN;

            % Compute mean particle number excluding NaNs
            smsAnnualDh(iTerm,iDh,:) = squeeze(mean(smsValues, 1, 'omitnan'));
        end % iDh
    end % iTerm
    
    % Process auxilliary terms
    % Define auxiliary terms that only have a single data layer
    singleLayerTerms = [config.aux.idxDiatBiomass, config.aux.idxFlagelBiomass, ...
                        config.aux.idxCoccoBiomass, config.aux.idxPicoBiomass, ...
                        config.aux.idxFreshDiatCellQuota, config.aux.idxFreshFlagelCellQuota, ...
                        config.aux.idxFreshCoccoCellQuota, config.aux.idxFreshPicoCellQuota, ...
                        config.aux.idxNightDvmUpperBound, config.aux.idxNightDvmLowerBound, ...
                        config.aux.idxDayDvmUpperBound, config.aux.idxDayDvmLowerBound, ...
                        config.aux.idxEuphoticDepth];
                
    for iTerm = 1:nAuxTerms
        if ismember(iTerm,singleLayerTerms)
            % Direct assignment for single-layer terms
            auxAnnualDh(iTerm,1,:) = squeeze(auxAnnualProfile(iTerm,1,:));
        else
            for iDh = 1:nTargetDepths
                % Apply precomputed mask
                currentMask = squeeze(validDepthMask(:,:,iDh));
                auxValues = squeeze(auxAnnualProfile(iTerm,:,:));
                auxValues(~currentMask) = NaN;

                % Compute mean, ignoring NaNs
                auxAnnualDh(iTerm,iDh,:) = squeeze(mean(auxValues, 1, 'omitnan'));
            end
        end
    end
    
end % extractKeyDepthHorizonDataInModelledSmsAndAuxTerms