function testProductionReadyMRGCmosaic

    % Retrieve the source RGC mosaic
    theSourceMidgetRGCMosaic = retrieveSourceRGCMosaic();


    % Instantiate a compute-ready optimized mRGCMosaic located
    % at (x,y) = (1,0.5), with width = 0.4 degs and height = 0.2 degs
    mySmallMRGCmosaic = mRGCMosaic(theSourceMidgetRGCMosaic, ...
        'eccentricityDegs', [1 0.5], ...
        'sizeDegs', [0.7 0.7], ...
        'name', 'my small off-center mRGC mosaic');

    theOI = mySmallMRGCmosaic.multiFocalRTVFopticsAtPosition(...
        mySmallMRGCmosaic.eccentricityDegs);
    

    % Create the letter A scene
    font = fontCreate('A', 'Georgia', 20, 96);
    display = 'LCD-Apple';
    theStimulusScene = sceneCreate('letter', font, display);
    theStimulusScene = sceneSet(theStimulusScene,'wangular',1.5);

    % Compute the retinal image
    theOI = oiCompute(theOI, theStimulusScene);

    % Compute the cone mosaic response to the retinal image placed at the center of the MRGC mosaic
    theConeMosaicResponse = mySmallMRGCmosaic.inputConeMosaic.compute(...
        theOI, 'opticalImagePositionDegs', mySmallMRGCmosaic.eccentricityDegs);

    % Visualize the cone mosaic response
    mySmallMRGCmosaic.inputConeMosaic.visualize(...
        'activation', theConeMosaicResponse);


    mySmallMRGCmosaic.visualize('component', 'RF centers');

    % Compute the MRGCmosaic response
    theConeMosaicResponseTemporalSupportSeconds = [0];
    theMRGCresponse = mySmallMRGCmosaic.compute(theConeMosaicResponse, ...
                             theConeMosaicResponseTemporalSupportSeconds);

    mySmallMRGCmosaic.visualize('activation', theMRGCresponse);
    
    pause;



    % Instantiate a compute-ready optimized mRGCMosaic with the same size
    % and position as the sourceMidgetRGCMosaic
    myLargeMRGCmosaic = mRGCMosaic(theSourceMidgetRGCMosaic, ...
        'name', 'my large mRGC mosaic');

    nTrials = 3;
    nTimePoints = 12;
    nCones = myLargeMRGCmosaic.inputConesNum;
    theConeMosaicResponse = 5+randn(nTrials, nTimePoints, nCones);
    theConeMosaicResponse(:,4,:) = 10;
    theConeMosaicResponseTemporalSupportSeconds = 4+(1:nTimePoints)*30/1000;

    [theMRGCresponses, theMRGCresponseTemporalSupportSeconds] = ...
        myMRGCmosaic.compute(theConeMosaicResponse, ...
                             theConeMosaicResponseTemporalSupportSeconds, ...
                             'timeResolutionSeconds', 10/1000);

    figure(1)
    for iTrial = 1:size(theMRGCresponses,1)
        subplot(1,3,iTrial)
        imagesc(theMRGCresponseTemporalSupportSeconds, 1:myMRGCmosaic.rgcsNum, (squeeze(theMRGCresponses(iTrial,:,:)))');
        xlabel('time (seconds');
        ylabel('RGC #');
    end

end


% ======== HELPER FUNCTIONS ========

function theSourceMidgetRGCMosaic = retrieveSourceRGCMosaic()
    dropboxDir = midgetRGCMosaicInspector.localDropboxPath;
    frozenMidgetRGCMosaicsDir = 'productionMidgetRGCMosaics/frozenMosaicsCenterSpecific';
    sourceMidgetRGCMosaicFileName = 'MRGCmosaic_Ecc_0.0_0.0_sizeDegs_3.0_3.0_H1cellIndex_1_Frozen.mat';
    sourceMidgetRGCMosaicFileName = fullfile(...
        dropboxDir,...
        frozenMidgetRGCMosaicsDir,...
        sourceMidgetRGCMosaicFileName);

    load(sourceMidgetRGCMosaicFileName, 'theMidgetRGCmosaic');
    theSourceMidgetRGCMosaic = theMidgetRGCmosaic;
    clear 'theMidgetRGCmosaic';
end
