function testProductionReadyMRGCmosaic

    dropboxDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';
    frozenMidgetRGCMosaicsDir = 'productionMidgetRGCMosaics/frozenMosaicsCenterSpecific';
    sourceMidgetRGCMosaicFileName = 'MRGCmosaic_Ecc_0.0_0.0_sizeDegs_3.0_3.0_H1cellIndex_1_Frozen.mat';
    sourceMidgetRGCMosaicFileName = fullfile(...
        dropboxDir,...
        frozenMidgetRGCMosaicsDir,...
        sourceMidgetRGCMosaicFileName);

    load(sourceMidgetRGCMosaicFileName, 'theMidgetRGCmosaic');

    % Instantiate a compute-ready optimized mRGCMosaic
    myMRGCmosaic = mRGCMosaic(theMidgetRGCmosaic, ...
        'name', 'my first RGC mosaic');

    nTrials = 3;
    nTimePoints = 12;
    nCones = myMRGCmosaic.inputConesNum;
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
