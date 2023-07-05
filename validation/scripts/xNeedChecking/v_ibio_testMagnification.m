function testMagnification

    s_initISET;
    
    % make a scene
    scene = sceneCreate('macbethd65');
    sceneAngularSizeInDeg = 1.0;
    sceneDistanceInMeters = 0.25;
    scene = sceneSet(scene,'wangular', sceneAngularSizeInDeg);
    scene = sceneSet(scene,'distance', sceneDistanceInMeters);
    
    % Computer optical image
    pupilDiameterInMillimeters = 4.0;  % 4 mm diameter
    pupilRadiusInMeters = pupilDiameterInMillimeters/2.0/1000.0;
    oi = oiCreate('human', pupilRadiusInMeters);
    
    oi = oiCompute(scene,oi);
    optics = oiGet(oi, 'optics');
    otfAll1 = opticsGet(optics,'otf data');
    otfSupport1 = oiGet(oi,'fsupport','cyclesPerDegree');
    size(otfSupport1)
    size(otfAll1)
    max(otfAll1(:))
    FullWidthHalfMaxBlurInMilliMeters = oiGet(oi, 'diffuser blur')*1000.0
    
    opticalRGBImage1     = oiGet(oi, 'rgb image');
    opticalImageWidthInMeters = oiGet(oi, 'width')
    
    
    pupilDiameterInMillimeters = 10.0;  % 4 mm diameter
    pupilRadiusInMeters = pupilDiameterInMillimeters/2.0/1000.0;
    oi = oiCreate('human', pupilRadiusInMeters);
    
    %waveSampling = 380:10:780;
    %oi  = oiSet(oi, 'wave', waveSampling);
    oi = oiCompute(scene,oi);
    optics = oiGet(oi, 'optics');
    otfAll2 = opticsGet(optics,'otf data');
    max(otfAll2(:))
    otfSupport2 = oiGet(oi,'fsupport','cyclesPerDegree');
    
    opticalRGBImage2     = oiGet(oi, 'rgb image');
    opticalImageWidthInMeters = oiGet(oi, 'width')
    
    FullWidthHalfMaxBlurInMilliMeters = oiGet(oi, 'diffuser blur')*1000.0
    
    figure(1);
    subplot(2,2,1);
    imshow(opticalRGBImage1)
    
    subplot(2,2,2);
    imshow(opticalRGBImage2)
    
    subplot(2,2,3);
    mesh(otfSupport1(:,:,1),otfSupport1(:,:,2),fftshift(abs(otfAll1)))
    
    subplot(2,2,4);
    mesh(otfSupport2(:,:,1),otfSupport2(:,:,2),fftshift(abs(otfAll2)))
    
    pause
    
    % Compute the optical image
    oi = oiCompute(scene,oi);
    optics = oiGet(oi, 'optics');

    magnification = opticsGet(optics, 'magnification', sceneDistanceInMeters)
    focalLengthInMeters = opticsGet(optics, 'focal length')
    computedMagnification = focalLengthInMeters/(focalLengthInMeters-sceneDistanceInMeters)

    sceneWidthInMeters  = sceneGet(scene, 'width');
    opticalImageWidthInMeters = oiGet(oi, 'width');
    
    imageMagnification = -opticalImageWidthInMeters/sceneWidthInMeters
    
    
    
    FullWidthHalfMaxBlurInMilliMeters = oiGet(oi, 'diffuser blur')*1000.0
    
    
end