function visualConeCharacteristicRadiusDegs = analyzeVisuallyProjectedConeAperture(...
                 anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
                 hFig, videoOBJ, pdfFileName)

    
    theCentroidArcMin = RetinaToVisualFieldTransformer.estimateGeometry(...
        thePSFData.supportX, thePSFData.supportY, thePSFData.data);


    % Generate anatomical cone aperture (a Gaussian with Rc)
    anatomicalConeCharacteristicRadiusArcMin = anatomicalConeCharacteristicRadiusDegs * 60;
    [Xarcmin, Yarcmin] = meshgrid(thePSFData.supportX, thePSFData.supportY);
    theAnatomicalConeApertureMap = exp(-(((Xarcmin-theCentroidArcMin(1)))/anatomicalConeCharacteristicRadiusArcMin).^2) .* ...
                                   exp(-(((Yarcmin-theCentroidArcMin(2)))/anatomicalConeCharacteristicRadiusArcMin).^2);

    % Convolve cone aperture with the PSF
    theVisuallyProjectedConeApertureMap = conv2(theAnatomicalConeApertureMap, thePSFData.data, 'same');
    theVisuallyProjectedConeApertureMap = theVisuallyProjectedConeApertureMap  / max(theVisuallyProjectedConeApertureMap(:));

    

    % Fit a 2D Gaussian to the visually projected cone aperture and extract
    % the characteristic radius of that Gaussian
    theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
        thePSFData.supportX, thePSFData.supportY, ...
        theVisuallyProjectedConeApertureMap);

    visualConeCharacteristicRadiusDegs = 1/60*sum(theFittedGaussian.characteristicRadii.^2,2)/sqrt(2);
    
    if (isempty(hFig))   
        return;
    end

    figure(2); clf;
    subplot(2,2,1);
    imagesc(thePSFData.supportX, thePSFData.supportY, theAnatomicalConeApertureMap);
    axis ('image');
    title('the anatomical cone aperture')
    subplot(2,2,2);
    imagesc(thePSFData.supportX, thePSFData.supportY, thePSFData.data);
    axis ('image');
    title('thePSF');

    subplot(2,2,3);
    imagesc(thePSFData.supportX, thePSFData.supportY, theVisuallyProjectedConeApertureMap);
    title('Visually projected cone aperture')
    axis ('image');
   
    subplot(2,2,4);
    imagesc(thePSFData.supportX, thePSFData.supportY, theFittedGaussian.ellipsoidMap);
    title('Fitted Ellipsoid');
    axis ('image');

   
    m = (size(theAnatomicalConeApertureMap,1)-1)/2+1;
    figure(3); clf;
    subplot(2,2,1);
    plot(thePSFData.supportX, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5);
    axis ('square');
    set(gca, 'XLim', 5*[-1 1], 'XTick', -5:1:5);
    grid on
    title('the anatomical cone aperture')

    subplot(2,2,2);
    plot(thePSFData.supportX, thePSFData.data(m,:), 'k-', 'LineWidth', 1.5);
    axis ('square');
    title('thePSF');
    set(gca, 'XLim', 5*[-1 1], 'XTick', -5:1:5);
    grid on
    subplot(2,2,3);
    plot(thePSFData.supportX, theAnatomicalConeApertureMap(m,:), 'k-', 'LineWidth', 1.5); hold on
    plot(thePSFData.supportX, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5);
    title('Visually & anatomical projected cone aperture')
    axis ('square');
    set(gca, 'XLim', 5*[-1 1], 'XTick', -5:1:5);
    grid on
    subplot(2,2,4);
    plot(thePSFData.supportX, theVisuallyProjectedConeApertureMap(m,:), 'r-', 'LineWidth', 1.5); hold on
    plot(thePSFData.supportX, theFittedGaussian.ellipsoidMap(m,:), 'b--', 'LineWidth', 1.5);
    plot(thePSFData.supportX, theFittedGaussian.ellipsoidMap(m,:)/max(theFittedGaussian.ellipsoidMap(:)), 'b:', 'LineWidth', 1.5);
    
    title('Fitted Ellipsoid');
    axis ('square');
    set(gca, 'XLim', 5*[-1 1], 'XTick', -5:1:5);
    grid on
end
