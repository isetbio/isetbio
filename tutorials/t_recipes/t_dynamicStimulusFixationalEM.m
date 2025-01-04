% Usage of computing with dynamic stimulus and fixational eye movements
%

% History:
%    01/03/25  NPC  ISETBIO Team, Copyright 2025 Wrote it.


%% Initialize
ieInit;

%% Stimulus params
% Duration: 0.5 second
stimulusDurationSeconds = 0.5;

% Frame duration: based on a 120 Hz refresh rate
% The time resolution of the fixationalEM and the
% integration time of the cone mosaic also get set to this value
frameDurationSeconds = 1/120;

% Field of view: 0.3 degs
fovDegs = .3;

%% Generate the OIsequence
theOISequence = generateOpticalImageSequence(stimulusDurationSeconds, frameDurationSeconds, fovDegs);

%% Generate the cone mosaic
theConeMosaic = cMosaic(...
    'sizeDegs', [0.5 0.5], ...      
    'eccentricityDegs', [0 0], ... 
    'integrationTime', frameDurationSeconds ...    
    );

%% Generate fixational eye movements
nTrials = 3;
theFixationalEMObj = generateFixationalEyeMovements(stimulusDurationSeconds, frameDurationSeconds, nTrials, theConeMosaic);

%% Compute mosaic responses
[theNeuralResponses, ~, ~, ~, temporalSupportSeconds] = ...
        theConeMosaic.compute(theOISequence, ...
        'withFixationalEyeMovements', true);

%% Visualize everything
visualizeConeExcitationsStimulusModulationAndFixationalEMs(...
    theConeMosaic, theNeuralResponses, temporalSupportSeconds, ...
    theOISequence, fovDegs, theFixationalEMObj);


%
% SUPPORT ROUTINES
%
function visualizeConeExcitationsStimulusModulationAndFixationalEMs(...
    theConeMosaic, theNeuralResponses, temporalSupportSeconds, ...
    theOISequence, fovDegs, theFixationalEMOb)

    stimulusIlluminanceSequence = [];
    for iTimePoint = 1:theOISequence.length
        stimulusIlluminanceSequence(iTimePoint,:,:) = oiGet(theOISequence.frameAtIndex(iTimePoint), 'illuminance');
    end
    illuminanceRange = [min(stimulusIlluminanceSequence(:)) max(stimulusIlluminanceSequence(:))];

    oiPixelWidthDegs = oiGet(theOISequence.frameAtIndex(iTimePoint), 'wangular resolution');
    oiWidthPixels = oiGet(theOISequence.frameAtIndex(iTimePoint), 'cols');
    oiSupport = (1:oiWidthPixels)*oiPixelWidthDegs;
    oiSupport = oiSupport - mean(oiSupport);

    hFig = figure(1);
    set(hFig, 'Position', [10 10 1800 650]);

    axConeMosaicExcitation = subplot(1,3,1);
    axStim = subplot(1,3,2);
    axEyeMovementPath = subplot(1,3,3);

    % Cone mosaic activation range
    activationRange = [min(theNeuralResponses(:)) max(theNeuralResponses(:))];

    % Visualize each frame of the stimulus/response/fixational EM
    nTrials = size(theNeuralResponses,1);
    timeSamplesNum = size(theNeuralResponses,2);

    

    for iTrial = 1:nTrials
    for iTimePoint = 1:timeSamplesNum
        theMosaicResponse = squeeze(theNeuralResponses(iTrial, iTimePoint,:));
        theConeMosaic.visualize('figureHandle', hFig,...
            'axesHandle', axConeMosaicExcitation, ...
            'activation', theMosaicResponse, ...
            'activationRange', activationRange, ...
            'plotTitleFontSize', 20, ...
            'plotTitle', sprintf('cone mosaic excitations\ntrial: %d, time: %2.1f msec', iTrial,temporalSupportSeconds(iTimePoint)*1e3));
    
            
        imagesc(axStim, oiSupport, oiSupport, (squeeze(stimulusIlluminanceSequence(iTimePoint,:,:))-illuminanceRange(1))/(illuminanceRange(2) - illuminanceRange(1)), [0 1]);
        axis(axStim, 'image');
        set(axStim, 'FontSize', 20, 'Color', [0 0 0]);
        set(axStim, 'XLim', theConeMosaic.sizeDegs(1)*0.5*[-1 1], 'YLim', theConeMosaic.sizeDegs(2)*0.5*[-1 1]);
        title(axStim, sprintf('retinal illuminance\n time: %2.1f msec', temporalSupportSeconds(iTimePoint)*1e3));
        colormap(gray);
    


        plot(axEyeMovementPath, ...
            theFixationalEMOb.emPosArcMin(iTrial, 1:iTimePoint,1)/60, ...
            theFixationalEMOb.emPosArcMin(iTrial, 1:iTimePoint,2)/60, '-', ...
            'LineWidth', 1.5, 'Color', 'r');
        axis(axEyeMovementPath, 'equal');
        set(axEyeMovementPath, 'XLim', theConeMosaic.sizeDegs(1)*0.5*[-1 1], 'YLim', theConeMosaic.sizeDegs(2)*0.5*[-1 1]);
        set(axEyeMovementPath, 'FontSize', 20, 'Color', [0 0 0]);
        title(axEyeMovementPath, sprintf('fixational eye movement path\ntrial: %d, time: %2.1f msec', iTrial,temporalSupportSeconds(iTimePoint)*1e3));
        drawnow;
         
    end % for iTimePoint
    end % for iTrial

end



function fixationalEMObj = generateFixationalEyeMovements(stimDurationSeconds, frameDurationSeconds, nTrials, theConeMosaic)
    % Initialize
    fixationalEMObj = fixationalEM;              % Instantiate a fixationalEM object
    fixationalEMObj.microSaccadeType = 'none';   % No microsaccades, just drift
    
    % Compute number of eye movements
    eyeMovementsPerTrial = stimDurationSeconds/frameDurationSeconds;

    % Generate the em sequence for the passed cone mosaic,
    % which results in a time step equal to the integration time of theConeMosaic
    fixationalEMObj.computeForCmosaic(...
        theConeMosaic, eyeMovementsPerTrial,...
        'nTrials' , nTrials);

    % Set the fixational eye movements into the cone mosaic
    theConeMosaic.emSetFixationalEMObj(fixationalEMObj);
end


function theOISequence = generateOpticalImageSequence(stimulusDurationSeconds, frameDurationSeconds, fovDegs)

    sparams = struct();
    
    vparams(2) = vernierP;
    vparams(2).name = 'offset';
    vparams(2).bgColor = 0;
    vparams(1) = vparams(2);
    vparams(1).barWidth = 0.1;
    vparams(1).bgColor = 0.5;
    vparams(1).name = 'uniform';
    sparams.fov = fovDegs;
    
    stimFramesNum = round(stimulusDurationSeconds/frameDurationSeconds);
    stimWeights = ieScale(fspecial('gaussian', [1, stimFramesNum], max([1 round(stimFramesNum/3)])), 0, 1);
    sampleTimes = (0:(numel(stimWeights)-1))*frameDurationSeconds;
    
    theOISequence = oisCreate('vernier', 'add', stimWeights, ...
       'testParameters', vparams, 'sceneParameters', sparams, ...
       'sampleTimes', sampleTimes);

end