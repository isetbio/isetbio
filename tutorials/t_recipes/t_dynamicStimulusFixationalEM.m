% Usage of computing with dynamic stimulus and fixational eye movements
%

% History:
%    01/03/25  NPC  ISETBIO Team, Copyright 2025 Wrote it.


%% Initialize
ieInit;

%% Stimulus params
% frame duration: 10 msec
frameDurationSeconds = 10/1000;
% FOV: 0.2 degs
fovDegs = 0.2;

%% Generate the OIsequence
theOISequence = generateOpticalImageSequence(frameDurationSeconds, fovDegs);

%% Generate the cone mosaic
theConeMosaic = cMosaic(...
    'sizeDegs', [0.5 0.5], ...      
    'eccentricityDegs', [0 0], ... 
    'integrationTime', frameDurationSeconds ...    
    );

%% Generate fixational eye movements
nTrials = 3;
theFixationalEMObj = generateFixationalEyeMovements(theConeMosaic, theOISequence.timeAxis, nTrials);

%% Compute mosaic responses
[theNeuralResponses, ~, ~, ~, temporalSupportSeconds] = ...
        theConeMosaic.compute(theOISequence, ...
        'withFixationalEyeMovements', true);

%% Visualize everything
visualizeEverything(theOISequence, fovDegs, theConeMosaic, theNeuralResponses, temporalSupportSeconds, theFixationalEMObj);


%
% SUPPORT ROUTINES
%
function visualizeEverything(theOISequence, fovDegs, theConeMosaic, theNeuralResponses, temporalSupportSeconds, theFixationalEMOb)

    stimulusIlluminanceSequence = [];
    for iTimePoint = 1:theOISequence.length
        stimulusIlluminanceSequence(iTimePoint,:,:) = oiGet(theOISequence.frameAtIndex(iTimePoint), 'illuminance');
    end
    illuminanceRange = [min(stimulusIlluminanceSequence(:)) max(stimulusIlluminanceSequence(:))];


    hFig = figure(1);
    set(hFig, 'Position', [10 10 1800 650]);

    axConeMosaicExcitation = subplot(1,3,1);
    axStim = subplot(1,3,2);
    axEyeMovementPath = subplot(1,3,3);

    % Cone mosaic activation range
    activationRange = [min(theNeuralResponses(:)) max(theNeuralResponses(:))];

    nTrials = size(theNeuralResponses,1);
    timeSamplesNum = size(theNeuralResponses,2);

    for iTrial = 1:nTrials
    for iTimePoint = 1:timeSamplesNum
        theMosaicResponse = squeeze(theNeuralResponses(iTrial, iTimePoint,:));
        theConeMosaic.visualize('figureHandle', hFig,...
            'axesHandle', axConeMosaicExcitation, ...
            'activation', theMosaicResponse, ...
            'activationRange', activationRange, ...
            'plotTitle', sprintf('cone mosaic excitations\ntrial: %d, t: %2.4f msec', iTrial,temporalSupportSeconds(iTimePoint)*1e3));
    
        xStimSupport = linspace(-fovDegs*0.5, fovDegs*0.5, size(stimulusIlluminanceSequence,3));
        imagesc(axStim, xStimSupport, xStimSupport, (squeeze(stimulusIlluminanceSequence(iTimePoint,:,:))-illuminanceRange(1))/(illuminanceRange(2) - illuminanceRange(1)), [0 1]);
        axis(axStim, 'image');
        set(axStim, 'FontSize', 20, 'Color', [0 0 0]);
        set(axStim, 'XLim', theConeMosaic.sizeDegs(1)*0.5*[-1 1], 'YLim', theConeMosaic.sizeDegs(2)*0.5*[-1 1]);
        title(axStim, 'stimulus');
        colormap(gray);
    

        cla(axEyeMovementPath);
        plot(axEyeMovementPath, ...
            theFixationalEMOb.emPosArcMin(iTrial, 1:iTimePoint,1)/60, ...
            theFixationalEMOb.emPosArcMin(iTrial, 1:iTimePoint,2)/60, 'r-', 'LineWidth', 1.5);
        axis(axEyeMovementPath, 'equal');
        set(axEyeMovementPath, 'XLim', theConeMosaic.sizeDegs(1)*0.5*[-1 1], 'YLim', theConeMosaic.sizeDegs(2)*0.5*[-1 1]);
        set(axEyeMovementPath, 'FontSize', 20, 'Color', [0 0 0]);
        title(axEyeMovementPath, 'fixational eye movements');
        drawnow;
         
    end
    end

end



function fixationalEMObj = generateFixationalEyeMovements(theConeMosaic, theTimeAxis, nTrials)
    % Initialize
    fixationalEMObj = fixationalEM;              % Instantiate a fixationalEM object
    fixationalEMObj.microSaccadeType = 'none';   % No microsaccades, just drift
    
    % Specify parameters
    stimDuration = theTimeAxis(end)-theTimeAxis(1);
    frameDurationSeconds = theTimeAxis(2)-theTimeAxis(1);
    eyeMovementsPerTrial = stimDuration/frameDurationSeconds + 1;

    % Generate the em sequence
    fixationalEMObj.computeForCmosaic(...
        theConeMosaic, eyeMovementsPerTrial,...
        'nTrials' , nTrials);

    % Set the fixational eye movements into the cone mosaic
    theConeMosaic.emSetFixationalEMObj(fixationalEMObj);
end


function theOISequence = generateOpticalImageSequence(frameDurationSeconds, fovDegs)

    sparams = struct();
    vParams = struct();
    
    vparams(2) = vernierP;
    vparams(2).name = 'offset';
    vparams(2).bgColor = 0;
    vparams(1) = vparams(2);
    vparams(1).barWidth = 0.1;
    vparams(1).bgColor = 0.5;
    vparams(1).name = 'uniform';
    sparams.fov = fovDegs;
    
    stimWeights = ieScale(fspecial('gaussian', [1, 50], 15), 0, 1);
    sampleTimes = (0:numel(stimWeights)-1)*frameDurationSeconds;
    
    theOISequence = oisCreate('vernier', 'add', stimWeights, ...
       'testParameters', vparams, 'sceneParameters', sparams, ...
       'sampleTimes', sampleTimes);

end