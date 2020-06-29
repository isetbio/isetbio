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
    
    xRange = 5/60*[-1 1];
    figure(1); clf;
    for wIndex = 1:numel(wavelengthSupport)
        subplot(3,round(numel(wavelengthSupport)/3)+1,wIndex);
        imagesc(thePSFsupportDegs, thePSFsupportDegs, squeeze(thePSFs(:,:,wIndex)));
        axis 'square'
        set(gca, 'XLim', xRange, 'YLim', xRange);
        title(sprintf('%2.0f nm (%2.1f, %2.1f) degs', wavelengthSupport(wIndex), eccX, eccY));
    end
    colormap(gray);
end
