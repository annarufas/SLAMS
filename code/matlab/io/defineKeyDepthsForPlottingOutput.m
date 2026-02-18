function [keyDepthLevelsMonthlyData,keyDepthLevelsAnnualData] =... 
    defineKeyDepthsForPlottingOutput(config,zeuMonthly,choiceTypeGridDomain)

    seqZeuMonthly = NaN(12,config.nLocs,2);
    seqZeuAnnual = NaN(config.nLocs,2);
    
    if size(zeuMonthly,1) == config.nLats % zeuMonthly is a lat x lon x 12 array

        for iLoc = 1:config.nLocs
            if (choiceTypeGridDomain == 1) % global
                iLat = config.iyBb(iLoc);
                iLon = config.ixBb(iLoc);
            elseif (choiceTypeGridDomain == 2) % local
                iLat = iLoc;
                iLon = iLoc;
            end 
            zeuLocal = zeuMonthly(iLat,iLon,:);
            seqZeuMonthly(:,iLoc,1) = zeuLocal;
            seqZeuMonthly(:,iLoc,2) = config.errorFractionZeu.*zeuLocal;
            seqZeuAnnual(iLoc,1) = mean(zeuLocal,'omitnan');
            seqZeuAnnual(iLoc,2) = mean(config.errorFractionZeu.*zeuLocal,'omitnan');
        end
    
    else % zeuMonthly is a sequential array already
        
        seqZeuMonthly = zeuMonthly;
        for i = 1:2
            seqZeuAnnual(:,1) = mean(zeuMonthly(:,:,1),1,'omitnan');
            seqZeuAnnual(:,2) = mean(zeuMonthly(:,:,2),1,'omitnan');
        end
        
    end

    % Define depths
    monthlyZeuUpp = seqZeuMonthly(:,:,1) - seqZeuMonthly(:,:,2);
    monthlyZeuLow = seqZeuMonthly(:,:,1) + seqZeuMonthly(:,:,2);
    
    annualZeuUpp = seqZeuAnnual(:,1) - seqZeuAnnual(:,2);
    annualZeuLow = seqZeuAnnual(:,1) + seqZeuAnnual(:,2);

    zmidmesoUpp = config.zmidmesoDef - config.errorZmeso; 
    zmidmesoLow = config.zmidmesoDef + config.errorZmeso; 

    zmesoUpp = config.zmesoDef - config.errorZmeso;
    zmesoLow = config.zmesoDef + config.errorZmeso; 

    zbathyUpp = config.zbathyDef - config.errorZbathy;
    zbathyLow = config.zbathyDef + config.errorZbathy; 
    
    % Define depths for monthly data
    keyDepthLevelsMonthlyData = NaN(4,12,config.nLocs,2); % nKeyLevels x nMonths x nLocs x 2 (upp, low)
    keyDepthLevelsMonthlyData(1,:,:,1) = monthlyZeuUpp;
    keyDepthLevelsMonthlyData(1,:,:,2) = monthlyZeuLow;
    keyDepthLevelsMonthlyData(2,:,:,1) = repmat(zmidmesoUpp,12,config.nLocs); 
    keyDepthLevelsMonthlyData(2,:,:,2) = repmat(zmidmesoLow,12,config.nLocs);
    keyDepthLevelsMonthlyData(3,:,:,1) = repmat(zmesoUpp,12,config.nLocs);
    keyDepthLevelsMonthlyData(3,:,:,2) = repmat(zmesoLow,12,config.nLocs); 
    keyDepthLevelsMonthlyData(4,:,:,1) = repmat(zbathyUpp,12,config.nLocs); 
    keyDepthLevelsMonthlyData(4,:,:,2) = repmat(zbathyLow,12,config.nLocs); 

    % Define depths block for annual data
    keyDepthLevelsAnnualData = NaN(4,config.nLocs,2); % nKeyLevels x nLocs x 2 (upp, low)
    keyDepthLevelsAnnualData(1,:,1) = annualZeuUpp;
    keyDepthLevelsAnnualData(1,:,2) = annualZeuLow;
    keyDepthLevelsAnnualData(2,:,1) = repmat(zmidmesoUpp,1,config.nLocs);
    keyDepthLevelsAnnualData(2,:,2) = repmat(zmidmesoLow,1,config.nLocs); 
    keyDepthLevelsAnnualData(3,:,1) = repmat(zmesoUpp,1,config.nLocs); 
    keyDepthLevelsAnnualData(3,:,2) = repmat(zmesoLow,1,config.nLocs); 
    keyDepthLevelsAnnualData(4,:,1) = repmat(zbathyUpp,1,config.nLocs); 
    keyDepthLevelsAnnualData(4,:,2) = repmat(zbathyLow,1,config.nLocs); 
    
end % defineKeyDepthsForPlottingOutput