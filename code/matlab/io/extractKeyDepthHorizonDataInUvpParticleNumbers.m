function [uvpPnumMonthlyDh,uvpPnumMonthlyDhDepths] =...
    extractKeyDepthHorizonDataInUvpParticleNumbers(config,uvpPnumMonthlyProfile,uvpDepths)

    % Cap modelledDepths where uvpDepths stop
    modelledDepths = config.availImagSysDeployDepths;
    [~,idxLastObsDepth] = min(abs(modelledDepths - uvpDepths(end)));
    modelledDepths = modelledDepths(1:idxLastObsDepth);
    nModelledDepths = length(modelledDepths);

    [nSizeClasses,nObsDepths,nMonths,nLocs] = size(uvpPnumMonthlyProfile);
    uvpPnumMonthlyDh = NaN(nSizeClasses,nModelledDepths,nMonths,nLocs);
    uvpPnumMonthlyDhDepths = NaN(nSizeClasses,nModelledDepths,nMonths,nLocs);

    for iDh = 1:nModelledDepths
        for iSc = 1:nSizeClasses
            for iMonth = 1:nMonths
                for iLoc = 1:nLocs

                    % Find the first and second closest depth to modelledDepths
                    [~,sortedIdx] = sort(abs(uvpDepths - modelledDepths(iDh)));
                    firstIdx = sortedIdx(1);
                    secondIdx = sortedIdx(2);

                    if firstIdx > secondIdx
                        firstIdx = sortedIdx(2);
                        secondIdx = sortedIdx(1);
                    end

                    if ~isnan(firstIdx) && ~isnan(secondIdx) && firstIdx ~= 0 && secondIdx ~= 0
                        % Calculate the average observed depth
                        uvpPnumMonthlyDhDepths(iSc,iDh,iMonth,iLoc) =...
                            (uvpDepths(firstIdx)+uvpDepths(secondIdx))/2;
    
                        % Calculate the average observed value
                        uvpPnumMonthlyDh(iSc,iDh,iMonth,iLoc) = mean(...
                            uvpPnumMonthlyProfile(iSc,firstIdx:secondIdx,iMonth,iLoc));
                    end

                end
            end
        end
    end

end % extractKeyDepthHorizonDataInUvpParticleNumbers