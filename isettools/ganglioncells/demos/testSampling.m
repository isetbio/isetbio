function testSampling
    eccDegs = [0 12];
    subjectID = 10;
    pupilDiamMM = 3.0;
    wavelengthsListToCompute = 450:25:700;
    
    % Compute microns per degree factor for the target eccentricity
    sizeDegs = 0.1;
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, sqrt(sum(eccDegs.^2)));
    micronsPerDegree = sizeMicrons(1)/sizeDegs(1);

    % Compute eccentricity varying optics based on the Polans data
    [theOI, thePSF, psfSupportX, psfSupportY] = PolansOptics.oiForSubjectAtEccentricity(eccDegs, subjectID, ...
        pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, ...
        'noLCA', ~true, ...
        'subtractCentralRefraction', true);
    
    psfRangeArcMin = 10;
    figure(14); clf;
    for k = 1:numel(wavelengthsListToCompute)
        subplot(3,4,k);
        imagesc(psfSupportX, psfSupportY, squeeze(thePSF(:,:,k))); hold on;
        plot(psfRangeArcMin*[-1 1], [0 0], 'g-');
        plot([0 0], psfRangeArcMin*[-1 1], 'g-'); hold off;
        axis 'equal';
        set(gca, 'XLim', psfRangeArcMin*[-1 1], 'YLim', psfRangeArcMin*[-1 1]);
        title(sprintf('%d nm', wavelengthsListToCompute(k)));
    end
    colormap(gray)
    
end

