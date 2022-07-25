function visualConeCharacteristicRadiusDegs = analyzeEffectOfPSFonConeAperture(...
     anatomicalConeCharacteristicRadiusDegs, thePSFData, ...
     hFig, videoOBJ, pdfFileName)

    [Xarcmin, Yarcmin] = meshgrid(thePSFData.supportX, thePSFData.supportY);

    [theCentroid, RcX, RcY, theRotationAngle] = ...
        RetinaToVisualFieldTransformer.estimateGeometry(thePSFData.supportX, thePSFData.supportY, thePSFData.data);

    % Generate anatomical cone aperture (a Gaussian with Rc)
    anatomicalConeAperture = exp(-(((Xarcmin-theCentroid(1))/60)/anatomicalConeCharacteristicRadiusDegs).^2) .* ...
                             exp(-(((Yarcmin-theCentroid(2))/60)/anatomicalConeCharacteristicRadiusDegs).^2);

    % Convolve cone aperture with the PSF
    theVisuallyProjectedConeAperture = conv2(thePSFData.data, anatomicalConeAperture, 'same');
    theVisuallyProjectedConeAperture = theVisuallyProjectedConeAperture  / max(theVisuallyProjectedConeAperture(:));

    % Fit a 2D Gaussian to the visually projected cone aperture and extract
    % the characteristic radius of that Gaussian
    [visualConeCharacteristicRadiusDegs, visualConeCharacteristicMinorMajorRadiiDegs, theVisuallyProjectedConeApertureFittedGaussian, XYcenter] = ...
        RetinaToVisualFieldTransformer.fitGaussianEllipsoid(thePSFData.supportX, thePSFData.supportY, ...
        theVisuallyProjectedConeAperture);

    if (isempty(hFig))   
        return;
    end

    XLims = XYcenter(1) + 7*[-1 1];  % XLims
    YLims = XYcenter(2) + 7*[-1 1];  % YLims
    plotAnalysis(anatomicalConeAperture, thePSFData,...
        theVisuallyProjectedConeAperture,...
        theVisuallyProjectedConeApertureFittedGaussian, ...
        XLims, YLims, videoOBJ, pdfFileName)
end

function plotAnalysis(anatomicalConeAperture, thePSFData, ...
    theVisuallyProjectedConeAperture, ...
    theVisuallyProjectedConeApertureFittedGaussian, ...
    XLims, YLims, videoOBJ, pdfFileName)
    
    xRange = XLims(2)-XLims(1);
    yRange = YLims(2)-YLims(1);
    xyRange = max([xRange yRange]);
    if (xyRange < 10)
        xTick = -100:1:100;
    elseif (xyRange < 30)
        xTick = -100:2:100;
    elseif (xyRange < 50)
        xTick = -100:5:100;
    elseif (xyRange < 100)
        xTick = -200:10:200;
    else
        xTick = -400:20:400;
    end

    hFig = figure(10); clf;
    set(hFig, 'Color', [1 1 1], 'Position', [10 10 1600 400]);
    
    
    % Plot here
    cLUT= brewermap(1024, 'blues');
    zLevels = 0.05:0.15:1.0;

    % The cone aperture
    ax = subplot(1,4,1);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, anatomicalConeAperture, zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,'cone aperture');
    
    ax = subplot(1,4,2);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, thePSFData.data/max(thePSFData.data(:)), zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,'point spread function');

    ax = subplot(1,4,3);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, theVisuallyProjectedConeAperture, zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,sprintf('visually projected cone aperture\n conv(coneAperture, PSF)'));
    
    ax = subplot(1,4,4);
    contourf(ax,thePSFData.supportX, thePSFData.supportY, theVisuallyProjectedConeApertureFittedGaussian, zLevels);
    set(ax, 'XLim', XLims, 'YLim', YLims, 'FontSize', 14, 'Color', squeeze(cLUT(1,:)), 'XTick', xTick, 'YTick', xTick);
    axis(ax,'xy'); axis(ax, 'square'); 
    grid(ax, 'on');
    xlabel(ax,'arc min');
    xtickangle(ax, 0);
    title(ax,'fitted Gaussian ellipsoid');
    
    colormap(cLUT);

    drawnow;
    NicePlot.exportFigToPDF(pdfFileName, hFig, 300);

    if (~isempty(videoOBJ))
        videoOBJ.writeVideo(getframe(hFig));
    end

end