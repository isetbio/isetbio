function visualizeStimResponseMovie(responseTimeAxis, stimulusTimeAxis, stimSpatialParams, theStimulusRGBsequence, theMosaicResponses, theConeMosaic)

    global LCONE_ID
    global MCONE_ID
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
               'rowsNum', 1, ...
               'colsNum', 3, ...
               'heightMargin',  0.03, ...
               'widthMargin',    0.02, ...
               'leftMargin',     0.05, ...
               'rightMargin',    0.00, ...
               'bottomMargin',   0.04, ...
               'topMargin',      0.02);
           
    hFig = figure(56); clf;
    set(hFig, 'Position', [10 10 2048 750], 'Color', [1 1 1]);
    videoOBJ = VideoWriter('stimAndResponse', 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 60;
    videoOBJ.Quality = 100;
    videoOBJ.open();
    ax = subplot('Position', subplotPosVectors(1,1).v);
    theConeMosaic.visualizeGrid('axesHandle', ax, ...
        'apertureShape', 'disks', ...
        'ticksInVisualDegs', true, 'backgroundColor', [1 1 1]);

    cm = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    LMconeIndices = find((cm.coneTypes == LCONE_ID) | (cm.coneTypes == MCONE_ID));

    mosaicResponseRange = [
        min(min(theMosaicResponses(LMconeIndices,:)))
        max(max(theMosaicResponses(LMconeIndices,:)))];
    
    stimXaxis = (1:stimSpatialParams.pixelsNum)/stimSpatialParams.pixelsNum;
    stimXaxisDegs = (stimXaxis-mean(stimXaxis))*stimSpatialParams.fovDegs;
    xLimsDegs = theConeMosaic.fov(1)*0.5*[-1 1];
    yLimsDegs = theConeMosaic.fov(2)*0.5*[-1 1];
    
    stimInterval = responseTimeAxis(2) - responseTimeAxis(1);
    stimFrames = size(theStimulusRGBsequence,1);
    for k = 1:size(theMosaicResponses,2)
        ax = subplot('Position', subplotPosVectors(1,2).v);
        [~,stimFrame] = min(abs(responseTimeAxis(k)-stimInterval-stimulusTimeAxis));
        stimFrame = mod(stimFrame-1, stimFrames)+1;
        image(ax,stimXaxisDegs, stimXaxisDegs, squeeze(theStimulusRGBsequence(stimFrame,:,:,:)));
        title(sprintf('stimulus frame: %d', stimFrame));
        xlabel(ax, ''); ylabel(ax, '');
        axis 'image';
        set(ax, 'FontSize', 16, 'XTick', [], 'YTick', [], 'XLim', xLimsDegs, 'YLim', yLimsDegs);
        

        ax = subplot('Position', subplotPosVectors(1,3).v);
        activation = squeeze(theMosaicResponses(:,k));
        theConeMosaic.renderActivationMap(ax, activation, 'signalRange', mosaicResponseRange, ...
            'mapType', 'modulated disks', ...
            'backgroundColor', [0 0 0]);
        title(sprintf('time: %2.1f ms', responseTimeAxis(k)*1000));
        set(ax, 'FontSize', 16, 'XTick', [], 'YTick', []);
        xlabel(ax, ''); ylabel(ax, '');
        axis 'image'
        
         drawnow;
         videoOBJ.writeVideo(getframe(hFig));
    end
    
    videoOBJ.close();
end
