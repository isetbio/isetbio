function visualizePSFs(theOI, eccX, eccY)
    optics = oiGet(theOI, 'optics');
    wavelengthSupport = opticsGet(optics, 'otfwave');
    thePSFs = opticsGet(optics,'psf data');
    
    % Extract support in arcmin
    psfSupportMicrons = opticsGet(optics,'psf support','um');
    if (isfield(optics, 'micronsPerDegree'))
        micronsPerDegree = optics.micronsPerDegree;
    else
        focalLengthMeters = opticsGet(optics, 'focalLength');
        focalLengthMicrons = focalLengthMeters * 1e6;
        micronsPerDegree = focalLengthMicrons * tand(1);
    end

    xGridDegs = psfSupportMicrons{1}/micronsPerDegree;
    thePSFsupportDegs = xGridDegs(1,:);
    
    
    
    hFig = figure(1); clf;
    set(hFig,  'Color', [1 1 1]);
    colsNum = 3;
    rowsNum = 2;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', rowsNum, ...
               'colsNum', colsNum, ...
               'heightMargin',  0.03, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.05, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.04, ...
               'topMargin',      0.02);
           
    maxPSF = max(thePSFs(:));
    sampledWavelengths = [450 500 550 600 650 700];
    for k = 1:numel(sampledWavelengths)
        [~,sampledWindices(k)] = min(abs(wavelengthSupport-sampledWavelengths(k)));
    end

    supportInMicrons = true;
    if (supportInMicrons)
        spatialSupport = psfSupportMicrons{1};
        spatialSupport = spatialSupport(1,:);
        xRange = 30*[-1 1];
        spaceLabel = 'microns';
        spatialTicks = -30:10:30;
    else
        spatialSupport = thePSFsupportDegs*60;
        xRange = 5*[-1 1];
        spatialTicks = -6:2:6;
        spaceLabel = 'arc min';
    end
    

    for k = 1:numel(sampledWindices)
        wIndex = sampledWindices(k);
        row = floor((k-1)/colsNum)+1;
        col = mod(k-1,colsNum)+1;
        subplot('Position', subplotPosVectors(row,col).v);
        thePSF = squeeze(thePSFs(:,:,wIndex));
        imagesc(spatialSupport, spatialSupport, thePSF);
        hold on
        contour(spatialSupport, spatialSupport, squeeze(thePSFs(:,:,wIndex)), [0.1:0.1:0.9]*maxPSF, 'black', 'LineWidth', 0.7);
        axis 'square'
        set(gca, 'XLim', xRange, 'YLim', xRange, 'CLim', [0 maxPSF], 'ZLim', [0 maxPSF]);
        set(gca, 'XTick', spatialTicks, 'YTick', spatialTicks);
        if (row ==rowsNum)&&(col == 1)
           title(sprintf('%2.0f nm (ecc = %2.1f degs)', wavelengthSupport(wIndex), sqrt(eccX^2+eccY^2)));
           xlabel(sprintf('retinal space (%s)', spaceLabel));
        else
            title(sprintf('%2.0f nm', wavelengthSupport(wIndex)));
            set(gca, 'XTickLabel', {}, 'YTickLabel', {})
        end
        grid 'on'
        colormap(brewermap(1024, 'greys'));
        drawnow;
    end
    
    NicePlot.exportFigToPDF('PSFs', hFig, 300);
end
