function thePSFData = generateOpticsPSFdataForVisualization(theOI, visualizedWavelength, micronsPerDegree)

    optics = oiGet(theOI, 'optics');
    waves = opticsGet(optics, 'wave');

    psfSupportMicrons = opticsGet(optics,'psf support','um');
    xGridDegs = psfSupportMicrons{1}/micronsPerDegree;
    yGridDegs = psfSupportMicrons{2}/micronsPerDegree;

    thePSFData.supportXdegs = xGridDegs(1,:);
    thePSFData.supportYdegs = yGridDegs(:,1);

    [~,idx] = min(abs(visualizedWavelength-waves));
    targetWavelength = waves(idx);
    thePSFData.data = opticsGet(optics,'psf data',targetWavelength );
end