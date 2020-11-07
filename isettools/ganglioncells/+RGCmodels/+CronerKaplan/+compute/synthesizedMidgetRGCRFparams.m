function RFparams = synthesizedMidgetRGCRFparams(rgcPositionsDegs, rgcCenterConeInputsNum)
% Generate params for center and surround subregion

    % Make sure that inputs have the right dimensionality
    assert(size(rgcPositionsDegs,1) == numel(rgcCenterConeInputsNum), ...
        sprintf('Dimensions mismatch.'));
    
    % Compute ecc radii in degrees. NOTE: Here we are using only the horizontal eccentricity component
    % because the deconvolution files are computed along the horizontal meridian only (for now)
    rgcEccRadiiDegsToUseForDeconvolution = abs(rgcPositionsDegs(:,1));
    rgcEccRadiiDegsToUseForDeconvolution = sqrt(sum(rgcPositionsDegs.^2,2));
    
    % Compute deconvolution model
    deconvolutionModel = RGCmodels.CronerKaplan.compute.deconvolutionModel();
    
    % Compute interpolation indices and weights from tabulated data in nearby eccentricities.
    [interpolationEccIndices, interpolationEccWeights] = computeInterpolationIndices(...
        deconvolutionModel.tabulatedEccRadiiDegs, rgcEccRadiiDegsToUseForDeconvolution, 'center', 'eccentricities');
    
    % Memory allocation
    rgcsNum = size(rgcPositionsDegs,1);
    centerVisualCharacteristicRadiiDegs = zeros(rgcsNum,1);
    centerVisualPeakSensitivityAttenuation = zeros(rgcsNum,1);
    surroundRetinalCharacteristicRadiiDegs = zeros(rgcsNum,1);
    surroundVisualPeakSensitivityAttenuation = zeros(rgcsNum,1);
    
    % Use the deconvolutionModel.center to determine:
    % (a) the center's VISUAL characteristic radius based on the number of input cones to this cells' RF center, and 
    % (b) the center's visual peak sensitivity attenuation factor
    
    for RGCindex = 1:rgcsNum
        % Retrieve the interpolation indices and interpolation weights for
        % weighing characteristic radii from the 2 closest eccentricities
        neighboringEccIndices = interpolationEccIndices(RGCindex,:);
        weightsOfNeighboringEccs = interpolationEccWeights(RGCindex,:);
        
        % Retrieve number of input cones to this RGC RF center
        inputConesNum = rgcCenterConeInputsNum(RGCindex);
        
        % Retrieve the corresponding characteristic radii for the 2 closest eccentricities
        characteristicRadiusDegsAtNeihboringEccs = zeros(1,numel(neighboringEccIndices));
        peakSensitivitiesAtNeihboringEccs = zeros(1,numel(neighboringEccIndices));
        for k = 1:numel(neighboringEccIndices)
            % Find the coneInputIndex corresponding to the # of cones in this RF center
            coneInputsIndex = find(deconvolutionModel.center.coneInputsNum(neighboringEccIndices(k),:) == inputConesNum);
            characteristicRadiusDegsAtNeihboringEccs(k) = deconvolutionModel.center.characteristicRadiusDegs(neighboringEccIndices(k),coneInputsIndex);
            % Retrieve the peak sensitivities (ratio of retinal-to-visual peak sensitivity) for the 2 nearby eccentricities
            peakSensitivitiesAtNeihboringEccs(k) = deconvolutionModel.center.peakSensitivity(neighboringEccIndices(k),coneInputsIndex);
        end
        
        % Compute the cell's VISUAL characteristic radius as the weighted mean of the 2 nearby characteristic radii
        centerVisualCharacteristicRadiiDegs(RGCindex) = sum(characteristicRadiusDegsAtNeihboringEccs .* weightsOfNeighboringEccs,2);
       
        % Compute the attenuation in VISUAL peak sensitivity (due to the low-pass filtering of the optics, 
        % the visual image of a cone has lower peak amplitude than the retinal image of that cone)
        centerVisualPeakSensitivityAttenuation(RGCindex) = 1 / sum(peakSensitivitiesAtNeihboringEccs .* weightsOfNeighboringEccs,2);
    end % RGCindex
        
    % Use the Croner&Kaplan model to compute the center VISUAL peak sensitivity from the center's VISUAL characteristic radius 
    centerVisualPeakSensitivities = RGCmodels.CronerKaplan.constants.centerPeakSensitivityFromCharacteristicRadiusDegsForPcells(centerVisualCharacteristicRadiiDegs);
    
    % To extract the RETINAL peak sensitivity, multiply the VISUAL peak sensitivity by the centerVisualPeakSensitivityAttenuation (estimated by the deconvolution model)
    centerRetinalPeakSensitivities = centerVisualPeakSensitivities .* centerVisualPeakSensitivityAttenuation;
    
    % Having determined the center VISUAL characteristic radius, compute the surround VISUAL characteristic radius
    surroundVisualCharacteristicRadiiDegs = RGCmodels.CronerKaplan.constants.surroundRadiusFromCenterRadiusDegsForPcells(centerVisualCharacteristicRadiiDegs);

    % Next, we use the Croner&Kaplan model surround-to-center integrated VISUAL sensitivity ratio to compute the surround VISUAL peak sensitivity
    surroundVisualPeakSensitivities = RGCmodels.CronerKaplan.constants.surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells(rgcEccRadiiDegsToUseForDeconvolution) .* ...
        centerVisualPeakSensitivities ./ ((surroundVisualCharacteristicRadiiDegs./centerVisualCharacteristicRadiiDegs).^2);
    
    % In the final step, we compute the surround RETINAL characteristic radius and peak sensitivity 
    % using the deconvolutionModel.surround
    
    for RGCindex = 1:rgcsNum
        % Retrieve the interpolation indices and interpolation weights for
        % weighing characteristic radii from the 2 closest eccentricities
        neighboringEccIndices = interpolationEccIndices(RGCindex,:);
        weightsOfNeighboringEccs = interpolationEccWeights(RGCindex,:);
        
        % Retrieve number of input cones to this RGC RF center
        inputConesNum = rgcCenterConeInputsNum(RGCindex);
        
        % We have deconvolution properties for a number of surround radii. Lets find the 2 closest surrounds
        retinalCharacteristicRadiusEstimates = zeros(1,2);
        visualGainSensitivities = zeros(1,2);
        
        for k = 1:numel(neighboringEccIndices)
            % Find the coneInputIndex corresponding to the # of cones in this RF center
            coneInputsIndex = find(deconvolutionModel.center.coneInputsNum(neighboringEccIndices(k),:) == inputConesNum);
            
            % determine interpolation weights for 2 closest surround characteristic radii
            nonNanIndices = find(~isnan(squeeze(deconvolutionModel.surround.characteristicRadiusDegs(neighboringEccIndices(k),coneInputsIndex,:))));
            if (isempty(nonNanIndices))
                error('We have no surround data for centers with %d cone inputs', inputConesNum);
            end
            
            % The surround radii that were examined in the deconvolution model
            examinedVisualCharacteristicRadii = squeeze(deconvolutionModel.surround.characteristicRadiusDegs(neighboringEccIndices(k),coneInputsIndex,nonNanIndices));
            
            % Determine interpolation indices and weights for the 2 closest surround sizes
            [interpolationSurroundRadiiIndices, interpolationSurroundRadiiWeights] = computeInterpolationIndices(...
                examinedVisualCharacteristicRadii, surroundVisualCharacteristicRadiiDegs(RGCindex), 'surround', 'radii');
            
            % Mean (across 2 nearby examined surround radii) retinal characteristic radius
            examinedRetinalCharacteristicRadii = ...
                squeeze(deconvolutionModel.surround.nominalSurroundRetinalCharacteristicRadii(neighboringEccIndices(k),coneInputsIndex,interpolationSurroundRadiiIndices));
            retinalCharacteristicRadiusEstimates(k) = sum(examinedRetinalCharacteristicRadii' .* interpolationSurroundRadiiWeights,2);
            
            % Mean (across 2 nearby examined surround radii) peak sensitities  (ratio of retinal-to-visual peak sensitivity)
            examinedVisualPeakSensitivities = squeeze(deconvolutionModel.surround.peakSensitivity(neighboringEccIndices(k),coneInputsIndex,interpolationSurroundRadiiIndices));
            visualGainSensitivities(k) = sum(examinedVisualPeakSensitivities' .* interpolationSurroundRadiiWeights,2);
        end
        
        % Weighted retinal characteristic radius (according to the 2 neighboring eccentricities) estimates
        surroundRetinalCharacteristicRadiiDegs(RGCindex) = sum(retinalCharacteristicRadiusEstimates .* interpolationEccWeights(RGCindex,:),2);
        
        % Weighted attenuation in peak sensitivity 
        surroundVisualPeakSensitivityAttenuation(RGCindex) = 1 / sum(visualGainSensitivities .* weightsOfNeighboringEccs,2); 
    end % RGCindex
    
    % To extract the RETINAL peak sensitivity, multiply the Croner&Kaplan estimate of the visual peak sensitivity
    % by the surround VisualPeakSensitivityAttenuation (estimated by the deconvolution model)
    surroundRetinalPeakSensitivities = surroundVisualPeakSensitivities .* surroundVisualPeakSensitivityAttenuation;
    
    % Assemble RFparams struct
    RFparams = struct(...
        'rfEccRadiusDegsForDeconvolutionParams', rgcEccRadiiDegsToUseForDeconvolution, ...               % ecc of RGCs within the target patch
        'rfCenterPositionDegs', rgcPositionsDegs, ...
        'visual', struct(...                                                        % VISUAL RF properties
            'centerCharacteristicRadiiDegs', centerVisualCharacteristicRadiiDegs, ...             
            'surroundCharacteristicRadiiDegs', surroundVisualCharacteristicRadiiDegs, ...         
            'centerPeakSensitivities', centerVisualPeakSensitivities,...            
            'surroundPeakSensitivities', surroundVisualPeakSensitivities...          
            ), ...
         'retinal', struct(...                                                      % RETINAL RF properties  
             'centerCharacteristicRadiiDegs', nan(size(surroundRetinalCharacteristicRadiiDegs)), ...  % there is no retinal center characteristic radius
             'surroundCharacteristicRadiiDegs', surroundRetinalCharacteristicRadiiDegs, ...          
             'centerPeakSensitivities', centerRetinalPeakSensitivities,... 
             'surroundPeakSensitivities', surroundRetinalPeakSensitivities... 
            )...
         );
     
end


function  [interpolationIndices, interpolationWeights] = computeInterpolationIndices(...
        tabulatedValues, targetValues, subregionName, domainName)
    
    % Determine linear interpolation indices and interpolation values
    targetsNum = numel(targetValues);
    interpolationWeights = zeros(targetsNum,2);
    interpolationIndices = zeros(targetsNum,2);
    
    for targetIndex = 1:targetsNum
        
        % Find the index of the tabulated values that is greater or equal to the target value
        targetValue = targetValues(targetIndex);
        idxPos = find(tabulatedValues >= targetValue);
        if (isempty(idxPos))
            fprintf(2,'Deconvolution data for %s do not extend up to %2.3f %s. Max value: %2.3f.\n', ...
                subregionName,  targetValues(targetIndex), domainName, max(tabulatedValues));
            idxPos = length(tabulatedValues);
        end
        
        % Find the index of the tabulated values that is greater or equal to the target value
        idxNeg = find(tabulatedValues <= targetValue);
        if (isempty(idxNeg))
            fprintf(2,'Deconvolution data for %s  do not extend down to %2.3f %s. Min value: %2.3f.\n', ...
                subregionName, targetValues(targetIndex), domainName, min(tabulatedValues));
            idxNeg = 1;
        end
        
        % Intepolation eccentricity indices for this cell
        interpolationIndices(targetIndex,:) = [idxNeg(end) idxPos(1)];
        
        % Compute interpolation weights based on the distance of the target
        % value to the 2 intepolation values
        interpolationValues = tabulatedValues(interpolationIndices(targetIndex,:));
        interpolationValueRange = abs(diff(interpolationValues));
        
        if (interpolationIndices(targetIndex,1) == interpolationIndices(targetIndex,2))
            interpolationWeights(targetIndex,:) = [0.5 0.5];
        else 
            eccDiffs = abs(interpolationValues - targetValue);
            % Interpolation weights
            interpolationWeights(targetIndex,:) = [eccDiffs(2) eccDiffs(1)]/interpolationValueRange;
        end
        
        
%         fprintf('interpolation: cell at %2.2f, below/above: %2.2f, %2.2f: weights: (%2.2f,%2.2f)  \n', ...
%             targetValue, ...
%             tabulatedValues(interpolationIndices(targetIndex,1)), ...
%             tabulatedValues(interpolationIndices(targetIndex,2)), ...
%             interpolationWeights(targetIndex,1), interpolationWeights(targetIndex,2));
    end

end
