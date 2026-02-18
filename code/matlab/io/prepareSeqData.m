function seqProcessed = prepareSeqData(targetOutput,relUncertainty,dataType)

    % General function to add fractional uncertainties to fluxes or 
    % particle numbers (particle attributes is an annual array for which
    % error calcualtion is not a concern) and format data into a format 
    % suitable for further processing.
    %
    % Inputs:
    %   targetOutput     - Structure containing model outputs
    %   relUncertainty   - Fractional uncertainties (vector for flux, scalar for particle numbers)
    %   dataType         - 'flux', 'particleNumbers', 'particleAttributes', 'sms' or 'aux
    %
    % Output:
    %   seqProcessed     - data with uncertainties added (for flux and particle numbers)

    switch dataType
        case 'flux'
        
            % Extract flux components
            tracers = struct(...
                'POC', targetOutput.monthlyOrgCarbonFlux + targetOutput.monthlyTepFlux, ...
                'PIC', targetOutput.monthlyCalciteFlux, ...
                'BSi', targetOutput.monthlyOpalFlux ...
            );

            tracerNames = fieldnames(tracers);
            nTracers = numel(tracerNames);

            % Initialise output array: (depth x months x locations x 2 x tracerType)
            seqProcessed = NaN([size(tracers.POC), 2, nTracers]); 

            % Assign values and compute uncertainties
            for iTracer = 1:nTracers
                seqProcessed(:,:,:,1,iTracer) = tracers.(tracerNames{iTracer}); % Flux values
                seqProcessed(:,:,:,2,iTracer) = relUncertainty(iTracer);        % Uncertainty (constant per tracer)
            end 

        case 'particleNumbers'
        
            % Initialise output array: (size classes x depths x months x locations x 2)
            seqProcessed = NaN([size(targetOutput.monthlyParticleNumInSc), 2]);

            % Assign values and uncertainties
            seqProcessed(:,:,:,:,1) = targetOutput.monthlyParticleNumInSc; % Particle numbers
            seqProcessed(:,:,:,:,2) = relUncertainty;                   % Single uncertainty value
        
        case 'particleAttributes' % (main size classes x depths x attributes x locations)
            seqProcessed = targetOutput.annualParticleAvgAttByMainSc;
            
        case 'sms' % (nSmsTerms x depths x locations)
            seqProcessed = targetOutput.annualSMS;
            
        case 'aux' % (nAuxTerms x depths x locations)
            seqProcessed = targetOutput.annualAux;

        otherwise
            error('Invalid dataType. Use ''flux'', ''particleNumbers'', ''particleAttributes'', ''sms'' or ''aux''.');
    end
      
end % prepareSeqData
