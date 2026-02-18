function [muLimFactorPft, muNetPft, chlToCarbonPft] = calculatePhytoGrowthRate(...
    nitrate, silicicAcid, phosphate, temperature, parz, q10, tempRef,... 
    muMax0degPft, kNO3Pft, kPO4Pft, kSiDiat,  chlToCarbonMax, alphaChlSpecificPft,...
    muThreshold)

    % Preallocate outputs
    nDepths = size(parz,1);
    muLimFactorPft = zeros(4,nDepths);
    muNetPft       = zeros(4,nDepths);
    chlToCarbonPft = zeros(4,nDepths);

    % Convert PAR from W m-2 to umol photons m-2 s-1
    parz = parz.*1e6/(3.9e-19.*6.02e23);

    for iPft = 1:4

        % Temperature-dependent growth rate adjustment
        tempFunc = q10(iPft).^((temperature-tempRef)./10);
        muAtCurrentTemp = muMax0degPft(iPft).*tempFunc; % d-1

        % Nutrient limitation factor
        if (iPft == 1) % diatoms require silicic acid
            nutLim = min([nitrate./(nitrate+kNO3Pft(iPft)) (phosphate./(phosphate+kPO4Pft(iPft))) (silicicAcid./(silicicAcid+kSiDiat))],[],'all','omitnan'); % values between 0-1  
            %nutLim = min([nitrate./(nitrate+kNO3) (phosphate./(phosphate+kPO4))]).*(silicicAcid./(silicicAcid+kSiDiat));  
        else
            %nutLim = nitrate./(nitrate+kNO3);  % values between 0-1
            nutLim = min([nitrate./(nitrate+kNO3Pft(iPft)) (phosphate./(phosphate+kPO4Pft(iPft)))],[],'all','omitnan');
        end 

        % Chlorophyll-to-carbon ratio and initial slope of the P-I curve
        chlToCarbonPft(iPft,:) = chlToCarbonMax ./ ...
            (1 + (chlToCarbonMax .* alphaChlSpecificPft(iPft) .* parz) ./ ...
            (2 .* muAtCurrentTemp ./ (24*3600))); % g Chl (g C)-1
        
        initSlopePhotoIrrCurve = chlToCarbonPft(iPft,:)' .* alphaChlSpecificPft(iPft); % m2 (umol photons)-1

        % Light limitation factor
        kIrr = (muAtCurrentTemp./(3600*24))./initSlopePhotoIrrCurve; % umol photons m-2 s-1
        irrLim = parz./(sqrt(kIrr.^2+parz.^2)); % equivalent to irrLim = V/muT, where V = alpha*PARz*muT/ (sqrt(muT^2 + (alpha*PARz)^2))

        % Combined limitation factor
        muLimFactorPft(iPft,:) = nutLim.*irrLim;

        % Net growth rate (apply threshold)
        muNet = muAtCurrentTemp.*nutLim.*irrLim; 
        muNetPft(iPft,:) = min(muNet,muThreshold);

    end % iPft

end % calculatePhytoGrowthRate