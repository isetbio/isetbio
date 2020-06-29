function theOIsequence = generateOISequenceForDriftingGrating(theSceneFrames, theOI, stimSpatialParams, stimTemporalParams)
    sceneFramesNum = numel(theSceneFrames);
    sinusoidalPeriodSeconds = 1.0/stimTemporalParams.temporalFrequencyHz;
    framesPerSinusoidalPeriod = 360/stimSpatialParams.deltaPhaseDegs;
    frameDurationSeconds = sinusoidalPeriodSeconds/framesPerSinusoidalPeriod;

    % Generate the oiSequence
    oiTimeAxis = 0:frameDurationSeconds:(stimTemporalParams.stimDurationSeconds-frameDurationSeconds);
    nFrames = numel(oiTimeAxis);
    
    oiList = cell(1, nFrames);
    for oiFrame = 1:nFrames
        % Compute the retinal image
        theOI = oiCompute(theOI, theSceneFrames{mod(oiFrame-1,sceneFramesNum)+1});
        oiList{oiFrame} = theOI;
    end
    
    % Generate the oiSequence with the list of OIs
    theOIsequence = oiArbitrarySequence(oiList, oiTimeAxis);
end
