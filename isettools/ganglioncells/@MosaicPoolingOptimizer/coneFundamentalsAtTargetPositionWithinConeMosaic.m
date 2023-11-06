function customConeFundamentals = coneFundamentalsAtTargetPositionWithinConeMosaic(...
    theConeMosaic, theOptics, targetRegionPositionDegs, targetRegionSizeDegs)

    % Find cone indices within the target region
    if (numel(targetRegionSizeDegs) == 2)
        widthDegs = targetRegionSizeDegs(1);
        heightDegs = targetRegionSizeDegs(2);
    else
        widthDegs = targetRegionSizeDegs(1);
        heightDegs = targetRegionSizeDegs(1);
    end

    theStimulusRegion = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', targetRegionPositionDegs, ...
            'width', widthDegs, ...
            'height', heightDegs , ...
            'rotation', 0.0...
        ));
    targetConeIndices = theStimulusRegion.indicesOfPointsInside(theConeMosaic.coneRFpositionsDegs);

    idx = find(theConeMosaic.coneTypes(targetConeIndices) == cMosaic.LCONE_ID);
    indicesOfLconesWithinTargetRegion  = targetConeIndices(idx);

    idx = find(theConeMosaic.coneTypes(targetConeIndices) == cMosaic.MCONE_ID);
    indicesOfMconesWithinTargetRegion  = targetConeIndices(idx);

    idx = find(theConeMosaic.coneTypes(targetConeIndices) == cMosaic.SCONE_ID);
    indicesOfSconesWithinTargetRegion  = targetConeIndices(idx);

    fprintf('Target region contains %d L-cones, %d M-cones and %d S-cones\n', ...
        numel(indicesOfLconesWithinTargetRegion), numel(indicesOfMconesWithinTargetRegion), numel(indicesOfSconesWithinTargetRegion));

    % Initialize
    customConeFundamentals.wavelengthSupport = theConeMosaic.wave;
    customConeFundamentals.quantalExcitationSpectra = zeros(numel(theConeMosaic.wave),3);
    
    % Compute the cone mosaic's response to a series of monochromatic, 
    % spatially-uniform images with constant power in photons/sec-nm
    pixelSize = 64;
    
    for iMonoChromaticBand = 1:numel(customConeFundamentals.wavelengthSupport)
        % Set up a spatially uniform dummy scene
        % (black body spectrum of 5000 degK).
        scene  = sceneCreate('uniform bb', pixelSize, 5000, customConeFundamentals.wavelengthSupport);
    
        % Set the scene FOV to 1 deg
        fovDegs = 1.0;
        scene = sceneSet(scene,'fov', fovDegs);

        % Get the photons 
        photons = sceneGet(scene,'photons');

        % Zero photons at all wavelengths 
        photons = zeros(size(photons));

        % Target photons only at the current monomchromatic wavelength
        photonsPerSrM2NMSec = 1e25;
        photons(:,:,iMonoChromaticBand) = photonsPerSrM2NMSec * ones(pixelSize,pixelSize);
    
        % Rewrite the scene with the desired monochromatic spectrum
        scene = sceneSet(scene,'photons',photons);

        % Compute the retinal image of the monochromatic scene
        theMonochromaticOpticalImage = oiCompute(scene, theOptics);

        % Compute the L-, M-, and S-cone exciations to the monochromatic scene
        theConeExcitations = theConeMosaic.compute(...
            theMonochromaticOpticalImage, ...
            'opticalImagePositionDegs', targetRegionPositionDegs, ...
            'nTrials', 1);

        % Compute the quantal L-cone excitations at this wavelength
        customConeFundamentals.quantalExcitationSpectra(iMonoChromaticBand, cMosaic.LCONE_ID) = ...
            mean(theConeExcitations(indicesOfLconesWithinTargetRegion)) / photonsPerSrM2NMSec;

        % Compute the quantal M-cone excitations at this wavelength
        customConeFundamentals.quantalExcitationSpectra(iMonoChromaticBand, cMosaic.MCONE_ID) = ...
            mean(theConeExcitations(indicesOfMconesWithinTargetRegion)) / photonsPerSrM2NMSec;

        % Compute the quantal S-cone excitations at this wavelength
        customConeFundamentals.quantalExcitationSpectra(iMonoChromaticBand, cMosaic.SCONE_ID) = ...
            mean(theConeExcitations(indicesOfSconesWithinTargetRegion)) / photonsPerSrM2NMSec;
    end % iMonochromaticBand 
end