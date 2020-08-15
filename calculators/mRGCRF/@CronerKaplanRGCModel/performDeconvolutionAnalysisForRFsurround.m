function deconvolutionStruct = performDeconvolutionAnalysisForRFsurround(obj, deconvolutionStruct, ...
    conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits)
    
    deconvolutionStruct.surround = performDeconvolutionAnalysisForSurroundSubregion(obj,deconvolutionStruct.center, ...
         conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits);

end

function deconvolutionStruct = performDeconvolutionAnalysisForSurroundSubregion(obj,centerDeconvolutionStruct, ...
         conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits)

    % Deconvolution data is a container
    deconvolutionStruct.data = containers.Map(); 
    keyNames = keys(centerDeconvolutionStruct.data);
    
    % Generate the cone aperture profile
    coneApertureProfile = obj.generateConeApertureProfileForDeconvolution(thePSFsupportDegs, coneAperturesDegs);
    
    % Generate cone mask to mask out pixels outside the cone apertures
    inputConeIndices = 1:size(conePosDegs,1);
    coneMask = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegs);
    coneMask(coneMask>0.0001) = 1;
        
    for poolingSchemeIndex = 1:numel(keyNames)
        poolingSchemeName = keyNames{poolingSchemeIndex};
        
        % Visual RF center characteristic radius for this RF center pooling scheme (1-30 cones)
        centerVisualCharacteristicRadius = 0.5*(centerDeconvolutionStruct.data(poolingSchemeName).minVisualSigma+ ...
             centerDeconvolutionStruct.data(poolingSchemeName).maxVisualSigma);

        % Most likely visual RF surround characteristic radius
        targetSurroundVisualCharacteristicRadius = CronerKaplanRGCModel.surroundRadiusFromCenterRadiusDegs(centerVisualCharacteristicRadius);
        
        % Examine a number of retinal RF surround characteristic radii that
        % vary from [0.5 - 1.25] x visualRF surround
        surroundRetinalCharacteristicRadii = targetSurroundVisualCharacteristicRadius * [0.5 0.75 1.0 1.25];
        
        for surroundSchemeIndex = 1:numel(surroundRetinalCharacteristicRadii)
            surroundRetinalCharacteristicRadius = surroundRetinalCharacteristicRadii(surroundSchemeIndex);
            
            % compute cone weights with a Gaussian falloff
            coneWeights = coneWeightsForSurroundCharacteristicRadius(surroundRetinalCharacteristicRadius, conePosDegs);
            
            % Compute the retinal cone image (input cones with Gaussian activation)
            retinalConeImage = generateRetinalConeImage(coneWeights, conePosDegs, coneApertureProfile, thePSFsupportDegs);
            
            % Compute the visual cone image = retinal cone image * PSF
            visualConeImage = conv2(retinalConeImage, thePSF, 'same');
            
             % Integrate within cone apertures
            [retinalConeActivations, visualConeActivations, ...
                withinConeAperturesRetinalConeImage, ...
                withinConeAperturesVisualConeImage , ...
                withinConeAperturesGrid] = obj.integrateConeImageWithinConeAperturesForDeconvolution(retinalConeImage, visualConeImage, ...
                                            coneMask, conePosDegs, coneAperturesDegs, thePSFsupportDegs);

            % Normalize activations with respect to the retinal domain
            visualConeActivations = visualConeActivations / max(retinalConeActivations(:));
            retinalConeActivations = retinalConeActivations / max(retinalConeActivations(:));
            
            % For the surround always use the circular Gaussian
            functionName = 'circular Gaussian';
            
            fitTheContinuousConeImage = true;
            if (fitTheContinuousConeImage)
                % Fit the continous coneImage not the integrated discrete
                % responses (but only including pixels that are within the
                % cone apertures)
                [rfSigmasRetinal, rfGainRetinal,  hiResRetinalConeActivationMap, hiResPSFsupportDegs] = ...
                    obj.fitActivationMapForDeconvolution(functionName,  withinConeAperturesRetinalConeImage, withinConeAperturesGrid, ...
                    coneAperturesDegs, thePSFsupportDegs);

                [rfSigmasVisual, rfGainVisual,  hiResVisualConeActivationMap] = ...
                    obj.fitActivationMapForDeconvolution(functionName, withinConeAperturesVisualConeImage, withinConeAperturesGrid, ...
                    coneAperturesDegs, thePSFsupportDegs);
            else
                % Fit the integraded discrete cone activation maps
                [rfSigmasRetinal, rfGainRetinal, hiResRetinalConeActivationMap, hiResPSFsupportDegs] = ...
                    obj.fitActivationMapForDeconvolution(functionName, retinalConeActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs);

                [rfSigmasVisual, rfGainVisual, hiResVisualConeActivationMap] = ...
                    obj.fitActivationMapForDeconvolution(functionName, visualConeActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs);
            end
             
            % Log data in
            surroundRetinalGains(surroundSchemeIndex) = mean(rfGainRetinal);
            surroundVisualGains(surroundSchemeIndex) =  mean(rfGainVisual);
            surroundRetinalCharacteristicRadii(surroundSchemeIndex) = mean(rfSigmasRetinal);
            surroundVisualCharacteristicRadii(surroundSchemeIndex) = mean(rfSigmasVisual);
            
            % Visualize fits
            if (visualizeFits)
                visualizedSpatialSupportRange = 0.2*[-1 1];
                obj.visualizeFitsForDeconvolution(surroundSchemeIndex, thePSFsupportDegs, thePSF, ...
                    conePosDegs, retinalConeImage, visualConeImage, ...
                    retinalConeActivations, visualConeActivations, ...
                    hiResPSFsupportDegs, hiResRetinalConeActivationMap, hiResVisualConeActivationMap, ...
                    rfSigmasRetinal, rfGainRetinal, ...
                    rfSigmasVisual, rfGainVisual, ...
                    visualizedSpatialSupportRange);
            end
            
        end  % surround scheme index
        
        % Log data in
        deconvolutionStruct.data(poolingSchemeName) = struct(...
            'retinalGains', surroundRetinalGains, ...
            'visualGains',  surroundVisualGains, ...
            'retinalSigmas', surroundRetinalCharacteristicRadii, ...
            'visualSigmas', surroundVisualCharacteristicRadii ...
        );
    
    end % poolingSchemeIndex
    
end

function retinalConeImage = generateRetinalConeImage(coneWeights, conePosDegs, coneApertureProfile, thePSFsupportDegs)
    % A 2D grid of delta functions placed at the centers of cones
    retinalConeImage = zeros(length(thePSFsupportDegs), length(thePSFsupportDegs));
    for iCone = 1:size(conePosDegs,1)
         coneXpos = conePosDegs(iCone,1);
         coneYpos = conePosDegs(iCone,2);
         [~,ix] = min(abs(coneXpos - thePSFsupportDegs));
         [~,iy] = min(abs(coneYpos - thePSFsupportDegs));
         if ( (coneXpos>=thePSFsupportDegs(1)) && (coneXpos<=thePSFsupportDegs(end)) && ...
              (coneYpos>=thePSFsupportDegs(1)) && (coneYpos<=thePSFsupportDegs(end)) )
               retinalConeImage(iy,ix) = coneWeights(iCone);
         end
    end
    retinalConeImage = conv2(retinalConeImage, coneApertureProfile, 'same');
end


function coneWeights = coneWeightsForSurroundCharacteristicRadius(surroundRetinalCharacteristicRadius, conePosDegs)
    coneEccRadius = sqrt(sum(conePosDegs.^2,2));
    coneWeights = exp(-(coneEccRadius/surroundRetinalCharacteristicRadius).^2);
end