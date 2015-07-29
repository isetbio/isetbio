function demoOIwindowIssue

    % please uncomment one of the following 5 conditions
    conditionIndex = 'FAIL_1';
    % conditionIndex = 'FAIL_2';
    % conditionIndex = 'SUCCESS_1';
    % conditionIndex = 'SUCCESS_2';
    % conditionIndex = 'SUCCESS_3';

    switch conditionIndex
        case 'FAIL_1'
            % conditions for failure
            generatePlots1 = false;
            generatePlots2 = true;
            ieInitFlag = true;

        case 'FAIL_2'
            % conditions for failure
            generatePlots1 = true;
            generatePlots2 = false;
            ieInitFlag = true;

        case 'SUCCESS_1'
            % conditions for success
            generatePlots1 = false;
            generatePlots2 = false;
            ieInitFlag = true;

        case 'SUCCESS_2'
            generatePlots1 = true;
            generatePlots2 = true;
            ieInitFlag = true;

        case 'SUCCESS_3'
            generatePlots1 = false;
            generatePlots2 = true;
            ieInitFlag = false;
    end


    isetbioIrradianceSpdPhotons1 = computeIsomerizations(generatePlots1,ieInitFlag);
    isetbioIrradianceSpdPhotons2 = computeIsomerizations(generatePlots2,ieInitFlag);

    hFig = figure(10);
    set(hFig, 'Position', [10 10 990 230]);
    subplot(1,3,1);
    plot(isetbioIrradianceSpdPhotons1)
    title(sprintf('generatePlots = %d', generatePlots1));

    subplot(1,3,2);
    plot(isetbioIrradianceSpdPhotons2)
    title(sprintf('generatePlots = %d', generatePlots2))

    subplot(1,3,3);
    plot(isetbioIrradianceSpdPhotons1-isetbioIrradianceSpdPhotons2)
    title('difference')
end

function isetbioIrradianceSpdPhotons = computeIsomerizations(generatePlots, ieInitFlag)

    if (ieInitFlag), ieInit; end

   %% Create isetbio display 
    displayToTest = 'LCD-Apple';
    d = displayCreate(displayToTest);
    gammaTable = displayGet(d, 'gamma table');
    nInputLevels = size(gammaTable,1);

    RGBToTest = [0.3 0.73 0.42]';
    theRGBImage = ones(20,20,3);
    for i1 = 1:3
        theRGBImage(:,:,i1) = round((nInputLevels-1)*RGBToTest(i1));
    end
    sceneDegrees = 10;  % need large field
    scene = sceneFromFile(theRGBImage,'rgb',[],d);
    scene = sceneSet(scene,'fov', sceneDegrees);
    sceneSize = sceneGet(scene,'size');
    
    
    if (generatePlots)
        vcAddAndSelectObject(scene); sceneWindow;
    end

    %% Compute optial image, for delta function optics.
    oi = oiCreate('human');
    optics = oiGet(oi,'optics');
    optics = opticsSet(optics,'off axis method','skip');
    optics = opticsSet(optics,'otf method','skip otf');
    oi = oiSet(oi,'optics',optics);
    oi = oiCompute(scene,oi);
    oi = oiSet(oi,'fov',sceneDegrees);
    
    if (generatePlots)
        vcAddAndSelectObject(oi); oiWindow;
    end

    roiPixels = 10;
    rect = [sceneSize(2)/2,sceneSize(1)/2,roiPixels,roiPixels];
    roiRoiLocs = ieRoi2Locs(rect);
    isetbioIrradianceSpdPhotons  = oiGet(oi,'roi mean photons',  roiRoiLocs); 
    
end