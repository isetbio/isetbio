function [visualConeCharacteristicRadiusDegs, bestHorizontalResolutionRotationDegs] = analyzeVisuallyProjectedConeAperture(theConeMosaic, ...
                 anatomicalConeApertureDiameterDegs, thePSFData, simulateCronerKaplanEstimation, hFig)

    % Compute the centroid of the PSF
    theCentroidDegs = RetinaToVisualFieldTransformer.estimateGeometry(...
        thePSFData.psfSupportXdegs, thePSFData.psfSupportYdegs, thePSFData.vLambdaWeightedData);


    % Generate anatomical cone aperture (a Gaussian with Rc) map centered
    % on the centroid of the PSF
    maxRadius = anatomicalConeApertureDiameterDegs*30;
    dx = thePSFData.psfSupportXdegs(2)-thePSFData.psfSupportXdegs(1);
    spatialSupportXdegs = dx:dx:maxRadius;
    spatialSupportXdegs = [-fliplr(spatialSupportXdegs) 0 spatialSupportXdegs];
    spatialSupportYdegs = spatialSupportXdegs;


    spatialSupportDegs = [spatialSupportXdegs(:)-theCentroidDegs(1) spatialSupportYdegs(:)-theCentroidDegs(2)];
    theAnatomicalConeApertureMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeMosaic, ...
        [0 0], ...
        1.0, ...
        spatialSupportDegs);

    theAnatomicalConeApertureMap = theAnatomicalConeApertureMap / max(theAnatomicalConeApertureMap(:));
    % Convolve cone aperture with the PSF
    theVisuallyProjectedConeApertureMap = conv2(theAnatomicalConeApertureMap, thePSFData.vLambdaWeightedData, 'same');
    theVisuallyProjectedConeApertureMap = theVisuallyProjectedConeApertureMap  / max(theVisuallyProjectedConeApertureMap(:));


    bestHorizontalResolutionRotationDegs = [];
    if (simulateCronerKaplanEstimation)
        % Since RF parameters by Croner&Kaplan were based on gratings, to
        % approximate this (and to include features of this estimation) we sum
        % along the y-dimension of the visually projected cone aperture map
        % and subsequently fit a line-weighting function for a Gaussian 

        % Rotate theVisuallyProjectedConeApertureMap so as to maximize horizontal resolution
        [rotatedTargetVisualRFmap,bestHorizontalResolutionRotationDegs] = ...
            RetinaToVisualFieldTransformer.bestHorizontalResolutionRFmap(theVisuallyProjectedConeApertureMap, bestHorizontalResolutionRotationDegs);
        
        % Fit a 1D Gaussian line weighting function to the 1D profile 
        % (integration along the Y-dimension of the 2D visually projected
        % cone aperture map)
        theVisuallyProjectedConeApertureMapProfile = sum(rotatedTargetVisualRFmap,1);
        theFittedGaussianLineWeightingFunction = RetinaToVisualFieldTransformer.fitGaussianLineWeightingFunction(...
            spatialSupportXdegs, theVisuallyProjectedConeApertureMapProfile);

        % Return the characteristic radius in degrees
        visualConeCharacteristicRadiusDegs = theFittedGaussianLineWeightingFunction.characteristicRadius;
    
        if (isempty(hFig))   
            return;
        end

        figure(hFig); clf;
        theAnatomicalConeApertureProfile = sum(theAnatomicalConeApertureMap,1);
        plot(spatialSupportXdegs, theVisuallyProjectedConeApertureMapProfile, 'ks-'); hold on
        plot(spatialSupportXdegs,theFittedGaussianLineWeightingFunction.profile, 'rs-');
        plot(spatialSupportXdegs,theAnatomicalConeApertureProfile, 'bs-');

        xlabel('space (degs)');
        legend('data', 'fit', 'cone aperture');
        axis 'square'
        
    else
        % Fit a 2D Gaussian ellipsoid to the 2D visually projected cone
        % aperture map
        theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            spatialSupportXdegs, spatialSupportXdegs, ...
            theVisuallyProjectedConeApertureMap);
    
        % Return the characteristic radius in degrees
        visualConeCharacteristicRadiusDegs = sqrt(sum(theFittedGaussian.characteristicRadii.^2,2))/sqrt(2);
    

        if (isempty(hFig))   
            return;
        end

        spatialSupportXLims = [min(spatialSupportXdegs) max(spatialSupportXdegs)];
        spatialSupportYLims = [min(spatialSupportYdegs) max(spatialSupportYdegs)];
        figure(hFig); clf;
        subplot(2,2,1);
        imagesc(spatialSupportXdegs*60, spatialSupportYdegs*60, theAnatomicalConeApertureMap);
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims*60, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        title('the anatomical cone aperture');
    
        subplot(2,2,2);
        imagesc(thePSFData.psfSupportXdegs*60, thePSFData.psfSupportYdegs*60, thePSFData.thePSFData.vLambdaWeightedData);
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims*60, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        title('thePSF');
    
        subplot(2,2,3);
        imagesc(spatialSupportXdegs*60, spatialSupportYdegs*60, theVisuallyProjectedConeApertureMap);
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims*60, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        title('Visually projected cone aperture')
        
       
        subplot(2,2,4);
        imagesc(spatialSupportXdegs*60, spatialSupportXdegs*60, theFittedGaussian.ellipsoidMap);
        title('Fitted Ellipsoid');
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims*60, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        colormap(gray)
       
        m = (size(theAnatomicalConeApertureMap,1)-1)/2+1;

        figure(3); clf;
        subplot(2,2,1);
        plot(spatialSupportXdegs*60, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5);
        axis ('square');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on
        xlabel('arcmin');
        title('the anatomical cone aperture')
    
        subplot(2,2,2);
        plot(thePSFData.psfSupportXdegs*60, thePSFData.data(m,:), 'k-', 'LineWidth', 1.5);
        axis ('square');
        title('thePSF');
        xlabel('arcmin');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on

        subplot(2,2,3);
        plot(spatialSupportXdegs*60, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5); hold on
        plot(spatialSupportXdegs*60, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5);
        title('Visually & anatomical projected cone aperture')
        axis ('square');
        xlabel('arcmin');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on
    
        subplot(2,2,4);
        plot(spatialSupportXdegs*60, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5); hold on
        plot(spatialSupportXdegs*60, theFittedGaussian.ellipsoidMap(m,:), 'b--', 'LineWidth', 1.5);
        plot(spatialSupportXdegs*60, theFittedGaussian.ellipsoidMap(m,:)/max(theFittedGaussian.ellipsoidMap(:)), 'b:', 'LineWidth', 1.5);
        
        title('Fitted Ellipsoid');
        axis ('square');
        xlabel('arcmin');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on

        disp('here')
        pause
    end
end
