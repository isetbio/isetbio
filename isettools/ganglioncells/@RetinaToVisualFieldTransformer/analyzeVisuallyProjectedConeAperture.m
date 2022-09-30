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
        % along the y-dimension of the visually projected cone aperture map
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
        spatialSupportXLims = [min(thePSFData.supportXdegs) max(thePSFData.supportXdegs)];
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

        spatialSupportXLims = [min(thePSFData.supportXdegs) max(thePSFData.supportXdegs)];
        spatialSupportYLims = [min(thePSFData.supportYdegs) max(thePSFData.supportYdegs)];
        figure(hFig); clf;
        subplot(2,2,1);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, theAnatomicalConeApertureMap);
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        title('the anatomical cone aperture');
    
        subplot(2,2,2);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, thePSFData.data);
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        title('thePSF');
    
        subplot(2,2,3);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, theVisuallyProjectedConeApertureMap);
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        title('Visually projected cone aperture')
        
       
        subplot(2,2,4);
        imagesc(thePSFData.supportXdegs*60, thePSFData.supportYdegs*60, theFittedGaussian.ellipsoidMap);
        title('Fitted Ellipsoid');
        axis ('image');
        set(gca,'XLim', spatialSupportXLims*60, 'YLim', spatialSupportYLims, 'XTick', -5:1:5, 'YTick', -5:1:5);
        xlabel('arcmin');
        ylabel('arcmin');
        colormap(gray)
       
        m = (size(theAnatomicalConeApertureMap,1)-1)/2+1;

        figure(3); clf;
        subplot(2,2,1);
        plot(thePSFData.supportXdegs*60, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5);
        axis ('square');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on
        xlabel('arcmin');
        title('the anatomical cone aperture')
    
        subplot(2,2,2);
        plot(thePSFData.supportXdegs*60, thePSFData.data(m,:), 'k-', 'LineWidth', 1.5);
        axis ('square');
        title('thePSF');
        xlabel('arcmin');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on

        subplot(2,2,3);
        plot(thePSFData.supportXdegs*60, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5); hold on
        plot(thePSFData.supportXdegs*60, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5);
        title('Visually & anatomical projected cone aperture')
        axis ('square');
        xlabel('arcmin');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on
    
        subplot(2,2,4);
        plot(thePSFData.supportXdegs*60, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5); hold on
        plot(thePSFData.supportXdegs*60, theFittedGaussian.ellipsoidMap(m,:), 'b--', 'LineWidth', 1.5);
        plot(thePSFData.supportXdegs*60, theFittedGaussian.ellipsoidMap(m,:)/max(theFittedGaussian.ellipsoidMap(:)), 'b:', 'LineWidth', 1.5);
        
        title('Fitted Ellipsoid');
        axis ('square');
        xlabel('arcmin');
        set(gca, 'XLim', spatialSupportXLims*60, 'XTick', -5:1:5);
        grid on
    end
end
