function dStruct = estimateConeCharacteristicRadiusInVisualSpace(theConeMosaic, thePSFData, theTargetPositionDegs,...
    coneCharacteristicRadiusConversionFactor, simulateCronerKaplanEstimation)

    conesNum = numel(theConeMosaic.coneTypes);
    if (conesNum == 0)
        fprintf(2, 'The mosaic contain no cones at this eccentricity, skipping computation of cone aperture in visual space.\n')
    
        % Return struct
        dStruct.conesNumInRetinalPatch = 0;
        dStruct.anatomicalConeCharacteristicRadiusDegs = nan;
        dStruct.visualConeCharacteristicRadiusDegs = nan;
        return;
    end
    
    % Sort cones according to their distance to theTargetPosition
    coneDistancesFromTargetPosition = sqrt(sum(bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, theTargetPositionDegs).^2,2));
    [~,idx] = sort(coneDistancesFromTargetPosition, 'ascend');
    
    % Estimate mean anatomical cone aperture from the 6 closest (to the target position) cones
    conesNumToUse = min([conesNum 6]);
    meanConeApertureDegs = mean(theConeMosaic.coneApertureDiametersDegs(idx(1:conesNumToUse)));
    anatomicalConeCharacteristicRadiusDegs = coneCharacteristicRadiusConversionFactor * meanConeApertureDegs;

    hFig = figure(1); clf;
    visualConeCharacteristicRadiusDegs = analyzeVisuallyProjectedConeAperture(...
        anatomicalConeCharacteristicRadiusDegs, thePSFData, simulateCronerKaplanEstimation, hFig);
    
    % Return struct
    dStruct.conesNumInRetinalPatch = conesNum;
    dStruct.indicesOfConesSortedWithDistanceToTargetRFposition = idx;
    dStruct.anatomicalConeCharacteristicRadiusDegs = anatomicalConeCharacteristicRadiusDegs;
    dStruct.visualConeCharacteristicRadiusDegs = visualConeCharacteristicRadiusDegs;
end

function visualConeCharacteristicRadiusDegs = analyzeVisuallyProjectedConeAperture(...
                 anatomicalConeCharacteristicRadiusDegs, thePSFData, simulateCronerKaplanEstimation, hFig)

    % Compute the centroid of the PSF
    theCentroidDegs = RetinaToVisualFieldTransformer.estimateGeometry(...
        thePSFData.supportXdegs, thePSFData.supportYdegs, thePSFData.data);


    % Generate anatomical cone aperture (a Gaussian with Rc) map centered
    % on the centroid of the PSF
    [Xdegs, Ydegs] = meshgrid(thePSFData.supportXdegs, thePSFData.supportYdegs);
    theAnatomicalConeApertureMap = exp(-(((Xdegs-theCentroidDegs(1)))/anatomicalConeCharacteristicRadiusDegs).^2) .* ...
                                   exp(-(((Ydegs-theCentroidDegs(2)))/anatomicalConeCharacteristicRadiusDegs).^2);

    % Convolve cone aperture with the PSF
    theVisuallyProjectedConeApertureMap = conv2(theAnatomicalConeApertureMap, thePSFData.data, 'same');
    theVisuallyProjectedConeApertureMap = theVisuallyProjectedConeApertureMap  / max(theVisuallyProjectedConeApertureMap(:));

    
    

    if (simulateCronerKaplanEstimation)
        % Since RF parameters by Croner&Kaplan were based on gratings, to
        % approximate this (and to include features of this estimation) we sum
        % along the y-dimension of the visualluy projected cone aperture map
        % and subsequently fit a line-weighting function for a Gaussian 

        % Fit a 1D Gaussian line weighting function to the 1D profile 
        % (integration along the Y-dimension of the 2D visually projected
        % cone aperture map)
        theVisuallyProjectedConeApertureMapProfile = sum(theVisuallyProjectedConeApertureMap,1);
        theFittedGaussianLineWeightingFunction = RetinaToVisualFieldTransformer.fitGaussianLineWeightingFunction(...
            thePSFData.supportXdegs, theVisuallyProjectedConeApertureMapProfile);

        % Return the characteristic radius in degrees
        visualConeCharacteristicRadiusDegs = theFittedGaussianLineWeightingFunction.characteristicRadius;
    
        if (isempty(hFig))   
            return;
        end

        figure(hFig); clf;
        plot(thePSFData.supportXdegs, theVisuallyProjectedConeApertureMapProfile, 'k-'); hold on
        plot(thePSFData.supportXdegs,theFittedGaussianLineWeightingFunction.profile, 'r-');
        xlabel('space (degs)');
        legend('data', 'fit');
        axis 'square'
    else
        % Fit a 2D Gaussian ellipsoid to the 2D visually projected cone
        % aperture map
        theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            thePSFData.supportXdegs, thePSFData.supportYdegs, ...
            theVisuallyProjectedConeApertureMap);
    
        % Return the characteristic radius in degrees
        visualConeCharacteristicRadiusDegs = sqrt(sum(theFittedGaussian.characteristicRadii.^2,2))/sqrt(2);
    

        if (isempty(hFig))   
            return;
        end

        figure(hFig); clf;
        subplot(2,2,1);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, theAnatomicalConeApertureMap);
        axis ('image');
        set(gca,'XLim', 3*[-1 1], 'YLim', 3*[-1 1], 'XTick', -3:0.5:3, 'YTick', -3:0.5:3);
        xlabel('arcmin');
        ylabel('arcmin');
        title('the anatomical cone aperture');
    
        subplot(2,2,2);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, thePSFData.data);
        axis ('image');
        set(gca,'XLim', 3*[-1 1], 'YLim', 3*[-1 1], 'XTick', -3:0.5:3, 'YTick', -3:0.5:3);
        xlabel('arcmin');
        ylabel('arcmin');
        title('thePSF');
    
        subplot(2,2,3);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, theVisuallyProjectedConeApertureMap);
        axis ('image');
        set(gca,'XLim', 3*[-1 1], 'YLim', 3*[-1 1], 'XTick', -3:0.5:3, 'YTick', -3:0.5:3);
        xlabel('arcmin');
        ylabel('arcmin');
        title('Visually projected cone aperture')
        
       
        subplot(2,2,4);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, theFittedGaussian.ellipsoidMap);
        title('Fitted Ellipsoid');
        axis ('image');
        set(gca,'XLim', 3*[-1 1], 'YLim', 3*[-1 1], 'XTick', -3:0.5:3, 'YTick', -3:0.5:3);
        xlabel('arcmin');
        ylabel('arcmin');
        colormap(gray)
       
        m = (size(theAnatomicalConeApertureMap,1)-1)/2+1;
        figure(3); clf;
        subplot(2,2,1);
        plot(thePSFData.supportXdegs*60, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5);
        axis ('square');
        set(gca, 'XLim', 3*[-1 1], 'XTick', -5:0.5:5);
        grid on
        xlabel('arcmin');
        title('the anatomical cone aperture')
    
        subplot(2,2,2);
        plot(thePSFData.supportXdegs*60, thePSFData.data(m,:), 'k-', 'LineWidth', 1.5);
        axis ('square');
        title('thePSF');
        xlabel('arcmin');
        set(gca, 'XLim', 3*[-1 1], 'XTick', -5:0.5:5);
        grid on

        subplot(2,2,3);
        plot(thePSFData.supportXdegs*60, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5); hold on
        plot(thePSFData.supportXdegs*60, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5);
        title('Visually & anatomical projected cone aperture')
        axis ('square');
        xlabel('arcmin');
        set(gca, 'XLim', 3*[-1 1], 'XTick', -5:0.5:5);
        grid on
    
        subplot(2,2,4);
        plot(thePSFData.supportXdegs*60, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5); hold on
        plot(thePSFData.supportXdegs*60, theFittedGaussian.ellipsoidMap(m,:), 'b--', 'LineWidth', 1.5);
        plot(thePSFData.supportXdegs*60, theFittedGaussian.ellipsoidMap(m,:)/max(theFittedGaussian.ellipsoidMap(:)), 'b:', 'LineWidth', 1.5);
        
        title('Fitted Ellipsoid');
        axis ('square');
        xlabel('arcmin');
        set(gca, 'XLim', 3*[-1 1], 'XTick', -5:0.5:5);
        grid on
    end
end
