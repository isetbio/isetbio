function testEccVaryingOpticsAndConeMosaic
    eccDegs = [10 0];
    subjectID = 10;
    pupilDiamMM = 3.0;
    
    [theConeMosaic, thePSF, wavelengthSupport, spatialSupport] = eccVaryingOpticsAndConeMosaic(eccDegs, subjectID, pupilDiamMM);
    
    % Visualize PSF and cone mosaic
    cMap = brewermap(1024, '*blues');
    zLevels = logspace(log10(0.01), log10(1), 6);
    psfRangeArcMin = 4;
    
    figure(14); clf;
    for k = 1:numel(wavelengthSupport)
        theWavePSF = squeeze(thePSF(:,:,k)) / max(thePSF(:));
        ax = subplot(3,4,k);
        PolansOptics.renderPSF(ax, spatialSupport, spatialSupport, theWavePSF, ...
            psfRangeArcMin, zLevels, cMap, 0.8*[1 1 1], ...
            'superimposedConeMosaic', theConeMosaic);
        title(ax,sprintf('%d nm', wavelengthSupport(k)));
    end
    drawnow;
    
end

function [theConeMosaic, thePSF, wavelengthSupport, psfSupportX] = eccVaryingOpticsAndConeMosaic(eccDegs, subjectID, pupilDiamMM)
    
    wavelengthSupport = 450:25:700;
    
    % Compute microns per degree factor for the target eccentricity
    sizeDegs = 0.1;
    eccRadiusDegs = sqrt(sum(eccDegs.^2));
    sizeMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(sizeDegs, eccRadiusDegs);
    micronsPerDegree = sizeMicrons(1)/sizeDegs(1);

    % Compute spacing for cone mosaic appropriate for this eccentricity
    coneSpacingDegs = RGCmodels.Watson.compute.coneSpacingAlongMeridianInRightEyeVisualField(eccRadiusDegs, 'temporal meridian');
    coneSpacingMicrons = RGCmodels.Watson.convert.sizeVisualDegsToSizeRetinalMicrons(coneSpacingDegs, eccRadiusDegs);
    coneApertureMicrons = RGCmodels.Watson.constants.coneApertureToDiameterRatio * coneSpacingMicrons;
    
    % Generate a 0.2 x 0.2 deg regular hex cone mosaic
    theConeMosaic = coneMosaicHex(9, ...
        'fovDegs', 0.2, ...
        'micronsPerDegree', micronsPerDegree, ...
        'customLambda', coneSpacingMicrons, ...
        'customInnerSegmentDiameter', coneApertureMicrons ...
    );
    
    % Compute optics for eccentricity  based on the Polans data
    [theOI, thePSF, psfSupportX, psfSupportY] = PolansOptics.oiForSubjectAtEccentricity(eccDegs, subjectID, ...
        pupilDiamMM, wavelengthSupport, micronsPerDegree, ...
        'wavefrontSpatialSamples', 301, ...
        'noLCA', ~true, ...
        'subtractCentralRefraction', true);
    
end

