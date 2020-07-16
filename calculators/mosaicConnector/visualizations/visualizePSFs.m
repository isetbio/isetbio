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
    
    xRange = 5*[-1 1];
    
    hFig = figure(1); clf;
    set(hFig,  'Color', [1 1 1]);
    colsNum = 4;
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
    for wIndex = 1:2:numel(wavelengthSupport)
        row = floor((wIndex-1)/colsNum)+1;
        col = mod(wIndex-1,colsNum)+1;
        subplot('Position', subplotPosVectors(row,col).v);
        thePSF = squeeze(thePSFs(:,:,wIndex));
        imagesc(thePSFsupportDegs*60, thePSFsupportDegs*60, thePSF);
        hold on
        contour(thePSFsupportDegs*60, thePSFsupportDegs*60, squeeze(thePSFs(:,:,wIndex)), [0.1:0.1:0.9]*maxPSF, 'black', 'LineWidth', 1.0);
        axis 'square'
        set(gca, 'XLim', xRange, 'YLim', xRange, 'CLim', [0 maxPSF], 'ZLim', [0 maxPSF]);
        set(gca, 'XTick', [-6:2:6], 'YTick', [-6:2:6]);
        if (row ==rowsNum)&&(col == 1)
           title(sprintf('%2.0f nm (ecc = %2.1f degs)', wavelengthSupport(wIndex), sqrt(eccX^2+eccY^2)));
           xlabel('retinal space (arc min)');
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
