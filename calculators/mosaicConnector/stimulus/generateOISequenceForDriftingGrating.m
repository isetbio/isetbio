function theOIsequence = generateOISequenceForDriftingGrating(theSceneFrames, theOptics, stimTemporalParams)

    sinusoidalPeriodSeconds = 1.0/stimTemporalParams.temporalFrequencyHz;
    sceneFramesNum = numel(theSceneFrames);
    frameDurationSeconds = sinusoidalPeriodSeconds / sceneFramesNum;
    oiTimeAxis = 0:frameDurationSeconds:stimTemporalParams.stimDurationSeconds;
    nFrames = numel(oiTimeAxis);
    
    oiList = cell(1, nFrames);
    parfor oiFrame = 1:nFrames
        fprintf('Computing optical image for frame %d/%d\n', oiFrame, nFrames);
        % Compute the retinal image
        theOI = oiCompute(theSceneFrames{mod(oiFrame-1,sceneFramesNum)+1}, theOptics);
        oiList{oiFrame} = theOI;
    end
    
    % Generate the oiSequence with the list of OIs
    theOIsequence = oiArbitrarySequence(oiList, oiTimeAxis);
end
