function examineRelationshipBetweenRetinalImageAndPSForientation
    
    testSubjectRankOrder = 6;
    rankedSujectIDs = PolansOptics.constants.subjectRanking();
    testSubjectID = rankedSujectIDs(testSubjectRankOrder);
    subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

    targetEcc = [-6 -6];
    pupilDiamMM = 3.0;
    wavefrontSpatialSamples = 501;
    zeroCenterPSF = true;
    flipPSFUpsideDown = true;
    upSampleFactor = 1.0;
    noLCA = false;
    refractiveErrorDiopters = 0;

    c = cMosaic('sizeDegs', [0.5 0.5]);

    [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength] = ...
                PolansOptics.oiForSubjectAtEccentricity(testSubjectID, ...
                'right eye', targetEcc, pupilDiamMM, c.wave, c.micronsPerDegree, ...
                'wavefrontSpatialSamples', wavefrontSpatialSamples, ...
                'subtractCentralRefraction', subtractCentralRefraction, ...
                'zeroCenterPSF', zeroCenterPSF, ...
                'flipPSFUpsideDown', flipPSFUpsideDown, ...
                'upsampleFactor', upSampleFactor, ...
                'noLCA',noLCA, ...
                'refractiveErrorDiopters', refractiveErrorDiopters);

    idx = find(psfSupportWavelength == 550);
    thePSF = squeeze(thePSF(:,:,idx));
    hFig = figure(1);
    ax = subplot(2,2,1);
    imagesc(ax,psfSupportMinutesX/60, psfSupportMinutesY/60, thePSF/max(thePSF(:)));
    axis(ax, 'image');
    set(ax, 'XLim', 0.3*[-1 1], 'YLim', 0.3*[-1 1]);

    theTestScene = sceneCreate('point array', 32);
    theTestScene  = sceneSet(theTestScene, 'fov', max(psfSupportMinutesY(:))/60);
    ax = subplot(2,2,2);
    visualizeScene(sceneGet(theTestScene, 'rgbimage'), ...
        'spatialSupportInDegs', true, ...
        'crossHairsAtOrigin', true ...
        'axesHandle', ax, ...
        );
    xlabel(ax,sprintf('%2.2f degs', max(psfSupportMinutesY(:))/60));
    axis(ax, 'image');

    theOI = oiCompute(theTestScene, theOI);
    ax = subplot(2,2,3);
    visualizeOpticalImage(theOI,...
        'crossHairsAtOrigin', true, ...
        'axesHandle', ax, ...
        'displayRadianceMaps', false);

    theActivationMap = c.compute(theOI);
    ax = subplot(2,2,4);
    c.visualize('figureHandle', hFig, ...
        'axesHandle', ax, ...
        'activation', theActivationMap);

end