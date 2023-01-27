function testProductionReadyMRGCmosaic

    % Retrieve the source RGC mosaic
    theSourceMidgetRGCMosaic = retrieveSourceRGCMosaic();


    % Instantiate a compute-ready optimized mRGCMosaic located
    % at (x,y) = (1,0.5), with width = 0.4 degs and height = 0.2 degs
    mySmallMRGCmosaic = mRGCMosaic(theSourceMidgetRGCMosaic, ...
        'eccentricityDegs', [1 0.5], ...
        'sizeDegs', [0.6 0.6], ...
        'name', 'my small off-center mRGC mosaic');

    theOI = mySmallMRGCmosaic.multiFocalRTVFopticsAtPosition(...
        mySmallMRGCmosaic.eccentricityDegs);
    

    % Stimulus params
    sceneFOVdegs = [2 1];
    stimulusPixelsNum = 256;
    retinalImageResolutionDegs = max(sceneFOVdegs)/stimulusPixelsNum;
    

    % Generate display
    viewingDistanceMeters = 2.0;
    thePresentationDisplay = rfMappingStimulusGenerator.presentationDisplay(...
            mySmallMRGCmosaic.inputConeMosaic.wave, retinalImageResolutionDegs, ...
            viewingDistanceMeters);

    stimSizeDegs = max(sceneFOVdegs);
    pixelSizeDegs = retinalImageResolutionDegs;
    spatialSupportDegs = rfMappingStimulusGenerator.spatialSupport(...
        stimSizeDegs, pixelSizeDegs);
    

    contrastMask = rfMappingStimulusGenerator.contrastMask(...
        spatialSupportDegs, ...
        'zeroInCentralRegionWithSize', [0.15 0.15] ...
        );

    coneContrasts = [1 1 1];
    orientationsTested = 0;

    spatialFrequenciesTested = [0.5];

    % Where to position  the stimulus
    stimulusRetinalPositionDegs = mySmallMRGCmosaic.eccentricityDegs;

    % Stim params for the STF mapping
    stimParams = struct(...
            'backgroundLuminanceCdM2', 50.0, ...
            'backgroundChromaticity', [0.301 0.301], ...
            'coneContrasts', coneContrasts, ...
            'contrast', 0.75, ...
            'orientationsTested', orientationsTested, ...
            'spatialFrequenciesTested', spatialFrequenciesTested, ...
            'orientationDegs', 0, ...
            'spatialFrequencyCPD', spatialFrequenciesTested(1), ...
            'spatialPhaseIncrementDegs', 5, ...
            'pixelSizeDegs', pixelSizeDegs, ...
            'stimSizeDegs', stimSizeDegs, ...
            'retinalPositionDegs', stimulusRetinalPositionDegs, ...
            'contrastMask', contrastMask, ...
            'wavelengthSupport', displayGet(thePresentationDisplay, 'wave'), ...
            'viewingDistanceMeters', displayGet(thePresentationDisplay, 'viewing distance') ...
            );

    [stimSpatialModulationFrames, stimParams.spatialPhasesDegs] = ...
        rfMappingStimulusGenerator.driftingGratingFrames(stimParams);

    [theDriftingGratingFrameScenes, theNullStimulusScene] = ...
                rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                    thePresentationDisplay, stimParams, stimSpatialModulationFrames, ...
                    'validateScenes', false);

    hFig = figure(10); clf;
    set(hFig, 'Position', [10 10 1450 700]);
    ax = subplot(1,2,1);
    image(sceneGet(theNullStimulusScene, 'rgbimage'));
    axis 'image'
    title('background scene')

    ax = subplot(1,2,2);
    image(sceneGet(theDriftingGratingFrameScenes{1}, 'rgbimage'));
    axis 'image'
    title('test stimulus (frame 1)')

    % Compute the retinal image of the background
    theBackgroundOI = oiCompute(theOI, theNullStimulusScene);

    % Compute the cone mosaic response to the retinal image of the background stimulus
        theNullConeMosaicResponse = mySmallMRGCmosaic.inputConeMosaic.compute(...
            theBackgroundOI, 'opticalImagePositionDegs', stimParams.retinalPositionDegs);

    for frameIndex = 1:numel(theDriftingGratingFrameScenes)

        % Compute the retinal image of first stimulus frame
        theStimulusOI = oiCompute(theOI, theDriftingGratingFrameScenes{frameIndex});

        % Compute the cone mosaic response to the retinal image of the test stimulus placed at the center of the MRGC mosaic
        theTestConeMosaicResponse(1,frameIndex,:) = mySmallMRGCmosaic.inputConeMosaic.compute(...
            theStimulusOI, 'opticalImagePositionDegs', stimParams.retinalPositionDegs);

        % Compute the modulation in the cone mosaic response (contrast)
        theConeMosaicResponseModulation(1,frameIndex,:) = bsxfun(@times, ...
            bsxfun(@minus, theTestConeMosaicResponse(1,frameIndex,:), theNullConeMosaicResponse), ...
            1./theNullConeMosaicResponse);
    end

    % Compute the MRGCmosaic response
    theConeMosaicResponseTemporalSupportSeconds = ((1:numel(theDriftingGratingFrameScenes))-1)* 30/1000;
    [theMRGCresponses, theMRGCresponseTemporalSupportSeconds] = mySmallMRGCmosaic.compute(...
                                theConeMosaicResponseModulation, ...
                                theConeMosaicResponseTemporalSupportSeconds);



    theROI = regionOfInterest('shape', 'line', 'from', [1,0.25], 'to', [1,0.75], 'thickness', 0.01);
    samplingPoints = 4000;  % sample the perimeter of the ROI along 1000 points');
    pointsPerSample = 10;  % return up to 6 points for each sample along the perimeter');
    maxDistance = 0.4;     % points must be no further than 0.03 degs away from the closest perimeter sample');

    % Find cones around the ROI
    idx = theROI.indicesOfPointsAround(...
        mySmallMRGCmosaic.inputConeMosaic.coneRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
    maxConesNumForTemporalResponsePlotting = 20;
    skip = floor(numel(idx)/maxConesNumForTemporalResponsePlotting);
    idx = idx(1:skip:numel(idx));
    identifiedConeIndices = idx;


    % Find RGCs around the ROI
    idx = theROI.indicesOfPointsAround(...
        mySmallMRGCmosaic.rgcRFpositionsDegs, pointsPerSample, samplingPoints, maxDistance);
    maxRGCsNumForTemporalResponsePlotting = 20;
    skip = floor(numel(idx)/maxRGCsNumForTemporalResponsePlotting);
    idx = idx(1:skip:numel(idx));
    identifiedRGCIndices = idx;


    s = abs(squeeze(theConeMosaicResponseModulation(1,:,identifiedConeIndices)));
    lineConeResponseRange = max(s(:))*[-1 1];

    generateConeMosaicVideo = false;
    if (generateConeMosaicVideo)
        videoOBJ = VideoWriter('ConeResponseDriftingGrating', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
     
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1950 950]);
    
        ax1 = subplot('Position', [0.07 0.07 0.4 0.9]);
        ax2 = subplot('Position', [0.55 0.07 0.4 0.9]);
        
        for timeBin = 1:size(theConeMosaicResponseModulation,2)
            % Visualize the cone mosaic response
            mySmallMRGCmosaic.inputConeMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax1, ...
                'activation', theConeMosaicResponseModulation(1,timeBin,:), ...
                'activationRange', max(abs(theConeMosaicResponseModulation(:))) * [-1 1], ...
                'backgroundColor', [0 0 0], ...
                'labelConesWithIndices', identifiedConeIndices, ...
                'verticalActivationColorBar', true, ...
                'plotTitle', sprintf('time: %2.0f msec', 1000.0*theConeMosaicResponseTemporalSupportSeconds(timeBin)));
    
            if (timeBin > 1)
                theIdentifiedConeResponses = squeeze(theConeMosaicResponseModulation(1,1:timeBin,identifiedConeIndices));
                plot(ax2, theConeMosaicResponseTemporalSupportSeconds(1:timeBin), ...
                          theIdentifiedConeResponses, 'r-', 'LineWidth', 1.5);
            end

            axis(ax2, 'square');
            set(ax2, 'YLim', lineConeResponseRange, 'XLim', [theConeMosaicResponseTemporalSupportSeconds(1) theConeMosaicResponseTemporalSupportSeconds(end)]);
            set(ax2, 'YTick', -1:0.1:1);
            set(ax2, 'FontSize', 16);
            xlabel(ax2, 'time (seconds)')
            ylabel(ax2, 'cone modulation');
            
    
            drawnow;
            videoOBJ.writeVideo(getframe(hFig));
        end
    
        videoOBJ.close()
    end



   
    % Visualize the mRGC mosaic response
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 2048 800]);
    
    ax1 = subplot('Position', [0.07 0.07 0.4 0.9]);
    ax2 = subplot('Position', [0.55 0.07 0.4 0.9]);

    mySmallMRGCmosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax1, ...
            'component', 'RF centers', ...
            'identifyInputCones', true);
        
    drawnow;


    hFig = figure(3); clf;
    set(hFig, 'Position', [10 10 2048 800], 'Color', [0 0 0]);
    ax1 = subplot('Position', [0.07 0.07 0.4 0.9]);
    ax2 = subplot('Position', [0.55 0.07 0.4 0.9]);
    
    activationColorMap = gray(512);
    activationColorMap(:,2) = 0;
    activationColorMap(:,3) = 0;
    redColorMap = activationColorMap;
    activationColorMap = gray(512);
    activationColorMap(:,1) = 0;
    activationColorMap(:,2) = 0;
    blueColorMap = activationColorMap;

    activationColorMap = cat(1, flipud(redColorMap), blueColorMap);

    activationColorMap = gray(1024);
    videoOBJ = VideoWriter('MRGCresponseDriftingGrating', 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();


    maxResponse = 0.8*max(abs(squeeze(theMRGCresponses(:))));
    s = abs(squeeze(theMRGCresponses(1,:,identifiedRGCIndices)));
    lineRGCResponseRange = max(s(:))*[-1 1];

    for timeBin = 1:size(theMRGCresponses,2)

        % Visualize the mRGC mosaic
%         mySmallMRGCmosaic.visualize(...
%             'figureHandle', hFig, ...
%             'axesHandle', ax1, ...
%             'activation', theMRGCresponse(1, timeBin,:), ...
%             'activationRange', maxResponse * [-1 1], ...
%             'verticalActivationColorBar', true, ...
%             'activationColorMap', activationColorMap, ...
%             'backgroundColor', activationColorMap(512,:), ...
%             'identifyInputCones', true, ...
%             'plotTitle', sprintf('time: %2.0f msec', ...
%                                   1000.0*theMRGCresponseTemporalSupportSeconds(timeBin)));


      
        mySmallMRGCmosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax1, ...
            'activation', theMRGCresponses(1, timeBin,:), ...
            'activationRange', maxResponse * [-1 1], ...
            'verticalActivationColorBar', true, ...
            'identifyInputCones', false, ...
            'identifyInputCones', true, ...
            'labelRGCsWithIndices', identifiedRGCIndices, ...
            'colorbarTickLabelColor', [0.8 0.8 0.8], ...
            'activationColorMap', activationColorMap, ...
            'backgroundColor', [0 0 0]);

        set(ax1, 'XColor', [0.7 0.7 0.7], 'YColor', [0.7 0.7 0.7]);

        if (timeBin> 1)
            cla(ax2);
            hold(ax2, 'on')

            theIdentifiedRGCresponses = squeeze(theMRGCresponses(1,1:timeBin,identifiedRGCIndices));

            for iRGC = 1:numel(identifiedRGCIndices)
                theRGCindex = identifiedRGCIndices(iRGC);
                switch (mySmallMRGCmosaic.centerSubregionMajorityConeTypes(theRGCindex))
                    case cMosaic.LCONE_ID
                        color = [1 0 0];
                    case cMosaic.MCONE_ID
                        color = [0 1 0];
                end

                plot(ax2, theMRGCresponseTemporalSupportSeconds(1:timeBin), ...
                      theIdentifiedRGCresponses(:, iRGC), '-', 'Color', color, 'LineWidth', 1);
            end
            
            plot(theMRGCresponseTemporalSupportSeconds(timeBin)*[1 1], lineRGCResponseRange, 'c-', 'LineWidth', 1.5)
            axis(ax2, 'square');
            set(ax2, 'YLim', lineRGCResponseRange, 'XLim', [theMRGCresponseTemporalSupportSeconds(1) theMRGCresponseTemporalSupportSeconds(end)]);
            set(ax2, 'Color', [0 0 0], 'YTick', -1:0.1:1);
            set(ax2, 'FontSize', 16, 'XColor', [0.7 0.7 0.7], 'YColor', [0.7 0.7 0.7]);
            xlabel(ax2, 'time (seconds)')
            ylabel(ax2, 'RGC modulation');
        end

        drawnow;
        videoOBJ.writeVideo(getframe(hFig));

    end

    videoOBJ.close();
    

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
