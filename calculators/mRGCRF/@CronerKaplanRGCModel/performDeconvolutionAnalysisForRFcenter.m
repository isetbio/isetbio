function deconvolutionStruct = performDeconvolutionAnalysisForRFcenter(obj, conesNumInRFcenterTested, ...
    conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits)
    
    deconvolutionStruct.center = performDeconvolutionAnalysisForCenterSubregion(obj,conesNumInRFcenterTested, ...
        conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits);
    
end


function deconvolutionStruct = performDeconvolutionAnalysisForCenterSubregion(obj,conesNumInRFcenterTested, ...
    conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits)

    % Center radius is about 6-7 times smaller than surround so only analyze central cones
    eccRadiiDegs = sqrt(sum(conePosDegs.^2,2));
    maxEccRadius = max(eccRadiiDegs);
    analyzedEccRadius = 0.3*maxEccRadius;
    coneIDsIncludedInCenterAnalysis = find(eccRadiiDegs <= analyzedEccRadius);
    fprintf('Only analyzing %d of the %d cones for the RF center size\n', numel(coneIDsIncludedInCenterAnalysis), numel(eccRadiiDegs));
    conePosDegs = conePosDegs(coneIDsIncludedInCenterAnalysis,:);
    coneAperturesDegs = coneAperturesDegs(coneIDsIncludedInCenterAnalysis);
    
    % Deconvolution data is a container
    deconvolutionStruct.data = containers.Map();  
    
    % Generate the cone aperture profile
    coneApertureProfile = obj.generateConeApertureProfileForDeconvolution(thePSFsupportDegs, coneAperturesDegs);
    
    % For some numbers of cone inputs, e.g. 2-5, we examine how the PSF interacts with several possible combinations of 
    % these inputs. Below we generate these possible input combinations
    rfCenterPoolingSchemes = generateCenterPoolingSchemes(conesNumInRFcenterTested, coneAperturesDegs, conePosDegs);
    
    inputConeIndices = 1:size(conePosDegs,1);
    coneMask = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegs);
    coneMask(coneMask>0.0001) = 1;
            
    for poolingSchemeIndex = 1:numel(rfCenterPoolingSchemes)
        
        poolingSchemeName = sprintf('%d-coneInput',rfCenterPoolingSchemes{poolingSchemeIndex}.coneInputsNum);
        
        inputConeIndicesForAllCombinations = rfCenterPoolingSchemes{poolingSchemeIndex}.inputConeIndices;
           
        for inputConeCombinationIndex = 1:size(inputConeIndicesForAllCombinations,1)
            fprintf('Analyzing %d of %d spatial input schemes for a %d-cone center subregion.\n',  ...
                inputConeCombinationIndex, size(inputConeIndicesForAllCombinations,1), rfCenterPoolingSchemes{poolingSchemeIndex}.coneInputsNum);
      
            % Retrieve the indices of cones feeding into the RF center
            inputConeIndices = inputConeIndicesForAllCombinations(inputConeCombinationIndex,:);
            
            % Compute the retinal cone image (all input cones lit-up)
            retinalConeImage = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegs);
            
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
            
            switch (numel(inputConeIndices))
                case 1
                    functionName = 'circular Gaussian';
                case 3
                    functionName = 'circular Gaussian';
                case 7
                    functionName = 'circular Gaussian';
                otherwise
                    functionName = 'elliptical Gaussian';
            end
            
            fitTheContinuousConeImage = true;
            if (fitTheContinuousConeImage)
                % Fit the continous coneImage not the integrated discrete
                % responses (but only including pixels that are within the
                % cone apertures)
                [rfSigmasRetinal(inputConeCombinationIndex,:), rfGainRetinal(inputConeCombinationIndex),  ...
                    hiResRetinalConeActivationMap, hiResPSFsupportDegs] = ...
                    obj.fitActivationMapForDeconvolution(functionName,  withinConeAperturesRetinalConeImage, withinConeAperturesGrid, ...
                    coneAperturesDegs, thePSFsupportDegs);

                [rfSigmasVisual(inputConeCombinationIndex,:), rfGainVisual(inputConeCombinationIndex),  ...
                    hiResVisualConeActivationMap] = ...
                    obj.fitActivationMapForDeconvolution('elliptical Gaussian', withinConeAperturesVisualConeImage, withinConeAperturesGrid, ...
                    coneAperturesDegs, thePSFsupportDegs);
            else
                % Fit the integraded discrete cone activation maps
                [rfSigmasRetinal(inputConeCombinationIndex,:), rfGainRetinal(inputConeCombinationIndex),  ...
                    hiResRetinalConeActivationMap, hiResPSFsupportDegs] = ...
                    obj.fitActivationMapForDeconvolution(functionName, retinalConeActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs);

                [rfSigmasVisual(inputConeCombinationIndex,:), rfGainVisual(inputConeCombinationIndex),  ...
                    hiResVisualConeActivationMap] = ...
                    obj.fitActivationMapForDeconvolution('elliptical Gaussian', visualConeActivations, conePosDegs, coneAperturesDegs, thePSFsupportDegs);
            end
            
            % Visualize fits
            if (visualizeFits)
                visualizedSpatialSupportRange = 0.05*[-1 1];
                obj.visualizeFitsForDeconvolution(inputConeCombinationIndex, thePSFsupportDegs, thePSF, ...
                    conePosDegs, retinalConeImage, visualConeImage, ...
                    retinalConeActivations, visualConeActivations, ...
                    hiResPSFsupportDegs, hiResRetinalConeActivationMap, hiResVisualConeActivationMap, ...
                    rfSigmasRetinal(inputConeCombinationIndex,:), rfGainRetinal(inputConeCombinationIndex), ...
                    rfSigmasVisual(inputConeCombinationIndex,:), rfGainVisual(inputConeCombinationIndex), ...
                    visualizedSpatialSupportRange);
            end
            
        end %inputConeCombination
        
        % Log data in
        deconvolutionStruct.data(poolingSchemeName) = struct(...
            'retinalGain', mean(rfGainRetinal), ...
            'visualGain',  mean(rfGainVisual), ...
            'minRetinalSigma', mean(min(rfSigmasRetinal,[],2)), ...
            'maxRetinalSigma', mean(max(rfSigmasRetinal,[],2)), ...
            'minVisualSigma', mean(min(rfSigmasVisual,[],2)), ...
            'maxVisualSigma', mean(max(rfSigmasVisual,[],2)) ...
        );
    end % poolingSchemeIndex      
end


function retinalConeImage = generateRetinalConeImage(inputConeIndices, conePosDegs, coneApertureProfile, thePSFsupportDegs)
    % A 2D grid of delta functions placed at the centers of cones
    retinalConeImage = zeros(length(thePSFsupportDegs), length(thePSFsupportDegs));
    for iCone = 1:numel(inputConeIndices)
         coneXpos = conePosDegs(inputConeIndices(iCone),1);
         coneYpos = conePosDegs(inputConeIndices(iCone),2);
         [~,ix] = min(abs(coneXpos - thePSFsupportDegs));
         [~,iy] = min(abs(coneYpos - thePSFsupportDegs));
         if ( (coneXpos>=thePSFsupportDegs(1)) && (coneXpos<=thePSFsupportDegs(end)) && ...
              (coneYpos>=thePSFsupportDegs(1)) && (coneYpos<=thePSFsupportDegs(end)) )
               retinalConeImage(iy,ix) = 1;
         end
    end
    retinalConeImage = conv2(retinalConeImage, coneApertureProfile, 'same');
end


function poolingSchemes = generateCenterPoolingSchemes(conesNumInRFcenterTested, coneAperturesDegs, conePosDegs)

    poolingSchemes = cell(1, numel(conesNumInRFcenterTested));
    
    coneCharacteristicRadiusDegs = 0.5*mean(coneAperturesDegs)/3;
    
    for schemeIndex = 1:numel(conesNumInRFcenterTested)
        coneInputsNum = conesNumInRFcenterTested(schemeIndex);
        switch (coneInputsNum)
            case 1
                scenariosNum = 1;
                startingAngle = 0; deltaAngle = 0;
                offsetRadius = 0;
            case 2
                scenariosNum = 6;
                startingAngle = 30; deltaAngle = 60;
                offsetRadius = coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            case 3
                scenariosNum = 6;
                startingAngle = 0; deltaAngle = 60;
                offsetRadius = 2/sqrt(3.0)*coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            case 4
                scenariosNum = 6;
                startingAngle = 0; deltaAngle = 50;
                offsetRadius = 2/sqrt(3.0)*coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            case 5
                scenariosNum = 2;
                startingAngle = 0; deltaAngle = 180;
                offsetRadius = coneCharacteristicRadiusDegs/WatsonRGCModel.coneApertureToDiameterRatio*3;
            otherwise
                scenariosNum = 1;
                startingAngle = 0; deltaAngle = 0;
                offsetRadius = 0;
        end
        
        poolingCenters = zeros(scenariosNum,2);
        inputConeIndices = zeros(scenariosNum,coneInputsNum);
        for iAngle = 1:scenariosNum
            poolingCenters(iAngle,1) = offsetRadius*cosd(startingAngle+iAngle*deltaAngle);
            poolingCenters(iAngle,2) = offsetRadius*sind(startingAngle+iAngle*deltaAngle);
            [~,inputConeIndices(iAngle,:)] = pdist2(conePosDegs, poolingCenters(iAngle,:), 'euclidean', 'smallest', coneInputsNum);
        end
        poolingSchemes{schemeIndex} = struct;
        poolingSchemes{schemeIndex}.coneInputsNum = coneInputsNum; 
        poolingSchemes{schemeIndex}.poolingCenters = poolingCenters;
        poolingSchemes{schemeIndex}.inputConeIndices = inputConeIndices;
    end
end