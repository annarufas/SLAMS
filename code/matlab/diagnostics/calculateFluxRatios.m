function [seqFluxRatMonthlyDh,seqFluxRatAnnualDh] = calculateFluxRatios(config,...
    seqFluxMonthlyDh,seqFluxAnnualDh)

    % Initialise output arrays
    seqFluxRatMonthlyDh = NaN([size(squeeze(seqFluxMonthlyDh(:,:,:,1,1))),2,2]); % (depth x months x locations x 2 x tracerType)
    seqFluxRatAnnualDh = NaN([size(squeeze(seqFluxAnnualDh(:,:,:,1,1))),2,2]);
    
    pocToPicMonthlyDh = NaN(size(squeeze(seqFluxMonthlyDh(:,:,:,1,1)))); 
    bSiToPicMonthlyDh = NaN(size(pocToPicMonthlyDh));
    pocToPicMonthlyDhErr = NaN(size(pocToPicMonthlyDh)); 
    bSiToPicMonthlyDhErr = NaN(size(pocToPicMonthlyDh));
    
    pocToPicAnnualDh = NaN(size(squeeze(seqFluxAnnualDh(:,:,1,1)))); 
    bSiToPicAnnualDh = NaN(size(pocToPicAnnualDh));
    pocToPicAnnualDhErr = NaN(size(pocToPicAnnualDh)); 
    bSiToPicAnnualDhErr = NaN(size(pocToPicAnnualDh));

    % Constants for unit conversion
    pocFactor = 1e3 * config.molarMassCarbon;
    calciteFactor = 1e3 * config.molarMassCaCO3;
    bSiFactor = 1e3 * config.molarMassBiogenicSilica;

    % Extract values and uncertainties
    pocFlux = squeeze(seqFluxMonthlyDh(:,:,:,1,1)); 
    picFlux = squeeze(seqFluxMonthlyDh(:,:,:,1,2));
    bSiFlux = squeeze(seqFluxMonthlyDh(:,:,:,1,3));
    pocFluxErr = squeeze(seqFluxMonthlyDh(:,:,:,2,1)); 
    picFluxErr = squeeze(seqFluxMonthlyDh(:,:,:,2,2)); 
    bSiFluxErr = squeeze(seqFluxMonthlyDh(:,:,:,2,3)); 

    % Calculate valid indices where picFlux is positive (PIC is in all
    % calculated ratios)
    validMask = (picFlux > 0);

    % Compute ratios and uncertainties where picFlux is positive
    pocToPicMonthlyDh(validMask) = (pocFlux(validMask) / pocFactor) ./...
        (picFlux(validMask) / calciteFactor);
    bSiToPicMonthlyDh(validMask) = (bSiFlux(validMask) / bSiFactor) ./...
        (picFlux(validMask) / calciteFactor);

    % Propagate uncertainties only at valid indices
    if any(validMask, 'all')
        pocToPicMonthlyDhErr(validMask) = calculateUncertaintyInQuotient(...
            (pocFlux(validMask) / pocFactor), (picFlux(validMask) / calciteFactor),...
            (pocFluxErr(validMask) / pocFactor), (picFluxErr(validMask) / calciteFactor));

        bSiToPicMonthlyDhErr(validMask) = calculateUncertaintyInQuotient(...
            (bSiFlux(validMask) / bSiFactor), (picFlux(validMask) / calciteFactor),...
            (bSiFluxErr(validMask) / bSiFactor), (picFluxErr(validMask) / calciteFactor));
    end

    % Add to output array
    seqFluxRatMonthlyDh(:,:,:,1,1) = pocToPicMonthlyDh;
    seqFluxRatMonthlyDh(:,:,:,1,2) = bSiToPicMonthlyDh;
    seqFluxRatMonthlyDh(:,:,:,2,1) = pocToPicMonthlyDhErr;
    seqFluxRatMonthlyDh(:,:,:,2,2) = bSiToPicMonthlyDhErr;

    % Compute annual averages while omitting NaNs
    pocToPicAnnualDh = squeeze(mean(pocToPicMonthlyDh, 2, 'omitnan'));
    bSiToPicAnnualDh = squeeze(mean(bSiToPicMonthlyDh, 2, 'omitnan'));
    
    % Propagate errors from monthly values to annual values
    for iDh = 1:3
        vals = squeeze(pocToPicMonthlyDh(iDh,:,:));
        errs = squeeze(pocToPicMonthlyDhErr(iDh,:,:));
        if any(~isnan(errs(:)))
            pocToPicAnnualDhErr(iDh,:,:) = calculateUncertaintyInAverage(vals,errs);
        end
    end
    for iDh = 1:3
        vals = squeeze(bSiToPicMonthlyDh(iDh,:,:));
        errs = squeeze(bSiToPicMonthlyDhErr(iDh,:,:));
        if any(~isnan(errs(:)))
            bSiToPicAnnualDhErr(iDh,:,:) = calculateUncertaintyInAverage(vals,errs);
        end
    end
    
    % Add to output array
    seqFluxRatAnnualDh(:,:,1,1) = pocToPicAnnualDh;
    seqFluxRatAnnualDh(:,:,1,2) = bSiToPicAnnualDh;
    seqFluxRatAnnualDh(:,:,2,1) = pocToPicAnnualDhErr;
    seqFluxRatAnnualDh(:,:,2,2) = bSiToPicAnnualDhErr;

end % calculateFluxRatios