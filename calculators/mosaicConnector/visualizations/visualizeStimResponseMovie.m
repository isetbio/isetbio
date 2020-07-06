function visualizeStimResponseMovie(responseTimeAxis, stimulusTimeAxis, theStimulusRGBsequence, theMosaicResponses, theConeMosaic)

    global LCONE_ID
    global MCONE_ID
    
    figure(56); clf;

    ax = subplot(1,2,1);
    theConeMosaic.visualizeGrid('axesHandle', ax);

    cm = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    LMconeIndices = find((cm.coneTypes == LCONE_ID) | (cm.coneTypes == MCONE_ID));

    mosaicResponseRange = max(max(abs(theMosaicResponses(LMconeIndices,:))))*[-1 1];
    
    stimInterval = responseTimeAxis(2) - responseTimeAxis(1);
    for k = 1:size(theMosaicResponses,2)
        ax = subplot(1,3,2);
        [~,stimFrame] = min(abs(responseTimeAxis(k)-stimInterval-stimulusTimeAxis));
        image(ax,squeeze(theStimulusRGBsequence(stimFrame,:,:,:)));
        title(sprintf('stimulus frame: %d', stimFrame));
        set(ax, 'FontSize', 12);
        axis 'image';

        ax = subplot(1,3,3);
        activation = squeeze(theMosaicResponses(:,k));
        theConeMosaic.renderActivationMap(ax, activation, 'signalRange', mosaicResponseRange);
        title(sprintf('bin: %d, time: %2.1f ms', k, responseTimeAxis(k)*1000));
        drawnow;
    end
end
