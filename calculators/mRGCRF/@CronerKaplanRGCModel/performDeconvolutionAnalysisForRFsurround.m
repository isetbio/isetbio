function deconvolutionStruct = performDeconvolutionAnalysisForRFsurround(obj, centerDeconvolutionStruct, ...
    examinedConesNumInRFCenter, conePosDegs, coneAperturesDegs, thePSF, thePSFsupportDegs, visualizeFits, exportFig, ...
            quadrantName, subjectID, patchEccRadiusDegs)
    
    % Interpolate PSF by a factor of 2
    upsampleFactor = 1;  % 3 for the center
    % Zero pad with a 0.3 degs margin on each size
    if (abs(patchEccRadiusDegs) > 15)
        paddingMarginDegs = 1.6;           % --> 2.0 degs padding: CHECK IF OK up to 20 cones @ 20 degs of ecc
    elseif (abs(patchEccRadiusDegs) > 10)
        paddingMarginDegs = 1.4;           % --> 1.4 degs padding: CHECK IF OK up to 20 cones @ 15 degs of ecc
    elseif (abs(patchEccRadiusDegs) > 5)
        paddingMarginDegs = 1.2;           % 1.2 degs padding CHECKED OK up to 20 cones @ 10 degs of ecc
    elseif (abs(patchEccRadiusDegs) > 2)
        paddingMarginDegs = 0.7;           % 0.7 degs padding CHECKED OK up to 20 cones @ 5 degs of ecc
    else
        paddingMarginDegs = 0.3;           % 0.3 degs padding CHECKED OK up to 20 cones @ 2 degs of ecc
    end
    
    [thePSFHR, thePSFsupportDegsHR] = ...
        CronerKaplanRGCModel.interpolatePSF(thePSF, thePSFsupportDegs, upsampleFactor, paddingMarginDegs);
    clear 'thePSFsupportDegs'
    
    % XY grid for 2D fitting
    [thePSFsupportDegsHRXgrid, thePSFsupportDegsHRYgrid] = meshgrid(thePSFsupportDegsHR,thePSFsupportDegsHR);
    
    % Generate the cone aperture profile
    [coneApertureProfile, ~] = obj.generateConeApertureProfileForDeconvolution(thePSFsupportDegsHR, coneAperturesDegs);
  
    % Ensure that the corresponding center deconvolution data were obtained with same parameters
    assert(patchEccRadiusDegs == centerDeconvolutionStruct.metaData.patchEccRadiusDegs, ...
        sprintf('patchEccRadiusDegs in surround (%2.2f) is different than in center (%2.2f)', ...
        patchEccRadiusDegs, centerDeconvolutionStruct.metaData.patchEccRadiusDegs));
    
    assert(strcmp(quadrantName, centerDeconvolutionStruct.metaData.quadrantName), ...
        sprintf('quadrantName in surround (''%s''), is different than in center (''%s'').', ...
        quadrantName, centerDeconvolutionStruct.metaData.quadrantName));
    
    assert(subjectID == centerDeconvolutionStruct.metaData.subjectID, ...
        sprintf('subjectID in surround (%d) is different than in center (%d).', ...
        subjectID, centerDeconvolutionStruct.metaData.subjectID));
    
    % DeconvolutionData is a container indexed by the number of cones in the RF center
    deconvolutionStruct.data = containers.Map();
    
    keyNames = keys(centerDeconvolutionStruct.data);
    for poolingSchemeIndex = 1:numel(examinedConesNumInRFCenter)
        % Get pooling scheme name
        poolingSchemeName = sprintf('%d-coneInput',examinedConesNumInRFCenter(poolingSchemeIndex));
        assert(ismember(poolingSchemeName, keyNames), ...
            sprintf('''%s'' is not a key in the centerDeconvolutionStruct.data container.', poolingSchemeName));
        
        % Visual RF center characteristic radius for this RF center pooling scheme (1-30 cones)
        centerVisualCharacteristicRadius = ...
            centerDeconvolutionStruct.data(poolingSchemeName).characteristicRadiusDegs;
        
        % Most likely visual RF surround characteristic radius
        targetSurroundVisualCharacteristicRadiusDegs = ...
            CronerKaplanRGCModel.surroundRadiusFromCenterRadiusDegs(centerVisualCharacteristicRadius);
        
        % Examine a number of retinal RF surround characteristic radii that
        % vary from [0.67 - 1.5] x visualRF surround
        surroundSizeVariations = (1:0.25:2.0);
        surroundSizeVariations = [1./fliplr(surroundSizeVariations) surroundSizeVariations(2:end)];
        
        if (visualizeFits)
            surroundSizeVariations = 1;
            fprintf(2, 'When visualizing surround deconvolution we only compute 1 surround size\n');
        end
        
        nominalSurroundRetinalCharacteristicRadiiDegs = targetSurroundVisualCharacteristicRadiusDegs * surroundSizeVariations;
        
        fittedPeakSensitivity = zeros(1, numel(nominalSurroundRetinalCharacteristicRadiiDegs));
        fittedCharacteristicRadiusDegs = zeros(1, numel(nominalSurroundRetinalCharacteristicRadiiDegs));
         
        parfor surroundSizeVariationIndex = 1:numel(nominalSurroundRetinalCharacteristicRadiiDegs)
            
            fprintf('Computing surround size variation %d/%d\n',surroundSizeVariationIndex, numel(nominalSurroundRetinalCharacteristicRadiiDegs));
            surroundRetinalCharacteristicRadius = nominalSurroundRetinalCharacteristicRadiiDegs(surroundSizeVariationIndex);
        
            % compute cone weights with a Gaussian falloff
            coneWeights = coneWeightsForSurroundCharacteristicRadius(surroundRetinalCharacteristicRadius, conePosDegs);
            
            % Compute the retinal cone image (input cones with Gaussian activation)
            retinalConeImage = generateRetinalConeImage(coneWeights, conePosDegs, coneApertureProfile, thePSFsupportDegsHR);
            
            % Compute the visual cone image = retinal cone image * PSF
            visualConeImage = conv2(retinalConeImage, thePSFHR, 'same');
            
            % Fit a circular Gaussian
            functionName = 'circular Gaussian';
            minCharacteristicRadiusDegs = centerVisualCharacteristicRadius;
            deltaX = thePSFsupportDegsHR(2)-thePSFsupportDegsHR(1);
            [fittedParams, rfFunction] = CronerKaplanRGCModel.fitElliptical2DGausianToRF(functionName, ...
                thePSFsupportDegsHRXgrid(:), thePSFsupportDegsHRYgrid(:), visualConeImage(:), ...
                deltaX, minCharacteristicRadiusDegs, [0 0]);

            % Compute fitted Gaussian
            xyData = zeros(numel(thePSFsupportDegsHRXgrid),2);
            xyData(:,1) = thePSFsupportDegsHRXgrid(:);
            xyData(:,2) = thePSFsupportDegsHRYgrid(:);
            N = numel(thePSFsupportDegsHR);
            fittedGaussian = reshape(rfFunction(fittedParams, xyData), [N N]);

            % Compute gain adjustment to match the areas of the
            % visualConeImage and the fittedGaussian
            visualSurroundArea = sum(visualConeImage(:));
            visualSurroundModelArea = sum(fittedGaussian(:));
            surroundGainAdjustmentToMatchAreas = visualSurroundArea/visualSurroundModelArea;

            % Adjust fittedGaussian and gain
            fittedGaussian = fittedGaussian * surroundGainAdjustmentToMatchAreas;
            fittedParams(1) = fittedParams(1) * surroundGainAdjustmentToMatchAreas;

            % Ensure that volumes match
            v1 = sum(visualConeImage(:));
            v2 = sum(fittedGaussian(:));
            assert(abs(v1-v2)/v1  <  1e-6, 'Volumes do not match');

            % Log data in
            fittedPeakSensitivity(surroundSizeVariationIndex) =  fittedParams(1);
            fittedCharacteristicRadiusDegs(surroundSizeVariationIndex) = fittedParams(4);

            if (visualizeFits)
                figNo = poolingSchemeIndex;
                visualizeAnalysis(figNo, subjectID, quadrantName, patchEccRadiusDegs, poolingSchemeName, ...
                    retinalConeImage, visualConeImage, fittedGaussian, thePSFsupportDegsHR, ...
                    surroundRetinalCharacteristicRadius, fittedCharacteristicRadiusDegs(surroundSizeVariationIndex), ...
                    exportFig);
            end
        end  % parfor surroundSizeVariationIndex
        fprintf('Finished with all surround variations for %s scheme (%d/%d) at ecc of %2.1f degs\n', poolingSchemeName, poolingSchemeIndex, numel(examinedConesNumInRFCenter), patchEccRadiusDegs);
        
        % Log data in
        deconvolutionStruct.data(poolingSchemeName) = struct(...
            'nominalSurroundRetinalCharacteristicRadii', nominalSurroundRetinalCharacteristicRadiiDegs, ...
            'characteristicRadiusDegs', fittedCharacteristicRadiusDegs, ...
            'peakSensitivity',  fittedPeakSensitivity ...
        );
    
    end % poolingSchemeIndex
end


function coneWeights = coneWeightsForSurroundCharacteristicRadius(surroundRetinalCharacteristicRadius, conePosDegs)
    coneEccRadius = sqrt(sum(conePosDegs.^2,2));
    coneWeights = exp(-(coneEccRadius/surroundRetinalCharacteristicRadius).^2);
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

function visualizeAnalysis(figNo, PolansSubjectID, quadrantName, patchEccRadiusDegs, poolingSchemeName, ...
            retinalConeImage, visualConeImage, fittedGaussian, thePSFsupportDegsHR, ...
            retinalCharacteristicRadius, visualCharacteristicRadius, exportFig)

    figWidthInches = 14;
    figHeightInches = 14;
    plotlabOBJ = setupPlotLab(0, figWidthInches, figHeightInches);
    
    hFig = figure(figNo); clf;
    set(hFig, 'Name', sprintf('subject %d, %s quadrant, ecc: %2.1f degs', ...
        PolansSubjectID, quadrantName, patchEccRadiusDegs));
    rowsNum = 2; colsNum = 2;
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum, ...
            'leftMargin', 0.05, ...
            'rightMargin', 0.00, ...
            'widthMargin', 0.08, ...
            'heightMargin', 0.09, ...
            'bottomMargin', 0.05, ...
            'topMargin', 0.01);
        
        
    zMax = max(retinalConeImage(:));
    
    ax = theAxesGrid{1,1}; 
    imagesc(ax,thePSFsupportDegsHR, thePSFsupportDegsHR, retinalConeImage);
    axis(ax, 'square');
    set(ax, 'CLim', [0 zMax], 'XLim', [-0.3 0.3], 'YLim', [-0.3 0.3]);
    title(ax, sprintf('retinal cone image (characteristic radius: %2.3f)', retinalCharacteristicRadius));
    colormap(ax,brewermap(1024, '*greys'));
    
    ax = theAxesGrid{1,2}; 
    imagesc(ax,thePSFsupportDegsHR, thePSFsupportDegsHR, visualConeImage);
    axis(ax, 'square');
    set(ax, 'CLim', [0 zMax], 'XLim', [-0.3 0.3], 'YLim', [-0.3 0.3]);
    title(ax,'visual cone image');
    colormap(ax,brewermap(1024, '*greys'));
    
    ax = theAxesGrid{2,2}; 
    imagesc(ax,thePSFsupportDegsHR, thePSFsupportDegsHR, fittedGaussian);
    axis(ax, 'square')
    set(ax, 'CLim', [0 zMax], 'XLim', [-0.3 0.3], 'YLim', [-0.3 0.3]);
    title(ax,sprintf('visual cone image (Gaussian model, characteristic radius: %2.3f)', visualCharacteristicRadius));
    colormap(ax,brewermap(1024, '*greys'));

    ax = theAxesGrid{2,1}; 
    midRow = floor(size(retinalConeImage,1)/2)+1;
    plot(ax,thePSFsupportDegsHR, retinalConeImage(midRow,:), 'r-', 'LineWidth', 1.5);
    hold(ax, 'on');
    plot(ax,thePSFsupportDegsHR, visualConeImage(midRow,:), 'b-', 'LineWidth', 1.5);
    plot(ax,thePSFsupportDegsHR, fittedGaussian(midRow,:), 'k-', 'LineWidth', 1.5);
    axis(ax,'square')
    set(ax, 'YLim', [0 zMax], 'XLim', [-0.3 0.3]);
    drawnow;
    
    % Export to PDF
    if (exportFig)
        isetbioPrefs = getpref('isetbio');
        isetbioRootDir = strrep(isetbioPrefs.validationRootDir, 'validation', '');
        exportDir = fullfile(isetbioRootDir, 'calculators/mosaicConnector/exports');
        pdfFileName = sprintf('Deconv_PolansSID_%d_%s_%2.1fdegs_%s_Surround',  ...
            PolansSubjectID, quadrantName, patchEccRadiusDegs, poolingSchemeName);
        plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, exportDir);
        %setupPlotLab(-1);
    end
    
end


function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 8, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'axesFontSize', 18, ...
                'axesFontAngle', 'italic', ...
                'axesLabelFontSizeMultiplier', 1.1, ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 


function OLD
        
        
    keyNames = keys(centerDeconvolutionStruct.data);
    
    % Deconvolution data is a container
    deconvolutionStruct.data = containers.Map(); 
    

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
        nominalSurroundRetinalCharacteristicRadii = targetSurroundVisualCharacteristicRadius * [0.5 0.75 1.0 1.25];
        
        for surroundSchemeIndex = 1:numel(nominalSurroundRetinalCharacteristicRadii)
            surroundRetinalCharacteristicRadius = nominalSurroundRetinalCharacteristicRadii(surroundSchemeIndex);
            
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
            'nominalSurroundRetinalCharacteristicRadii', nominalSurroundRetinalCharacteristicRadii, ...
            'retinalGains', surroundRetinalGains, ...
            'visualGains',  surroundVisualGains, ...
            'retinalSigmas', surroundRetinalCharacteristicRadii, ...
            'visualSigmas', surroundVisualCharacteristicRadii ...
        );
    
    end % poolingSchemeIndex
    
end



