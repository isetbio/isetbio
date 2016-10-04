function varargout = v_osTimeStep(varargin)
%
% Validate basic multi-optical image computations for different temporal dynamics. 
%
% This script checks multi-optical image computations with eye movements.
%
% NPC, ISETBIO TEAM, 2016
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;

% Assemble conditions to run.
condData = {};

conditionSet = 2;

if (conditionSet == 1)
    % Effects of varying the response time interval
    
    addCondition = true;
    if (addCondition)
    condData{numel(condData)+1} = struct(...
        'modulation', 0.1, ...                      % 10% photon modulation against background
        'stimulusSamplingInterval',  1/50, ...      % 50 Hz stimulus refresh
        'responseTimeInterval', 20/1000, ...        % 20 milliseconds
        'photonNoise', false, ...
        'osNoise', false);
    end

    addCondition = true;
    if (addCondition)
    condData{numel(condData)+1} = struct(...
        'modulation', 0.1, ...                      % 10% photon modulation against background
        'stimulusSamplingInterval',  1/50, ...      % 50 Hz stimulus refresh
        'responseTimeInterval', 5/1000, ...         % 5 milliseconds
        'photonNoise', false, ...
        'osNoise', false);
    end

    addCondition = true;
    if (addCondition)
    condData{numel(condData)+1} = struct(...
        'modulation', 0.1, ...                      % 10% photon modulation against background
        'stimulusSamplingInterval',  1/50, ...      % 50 Hz stimulus refresh
        'responseTimeInterval', 1/1000, ...         % 1 milliseconds
        'photonNoise', false, ...
        'osNoise', false);
    end
    
elseif (conditionSet == 2)
    % Effects of varying the noise and stimulus modulation
    addCondition = true;
    if (addCondition)
    condData{numel(condData)+1} = struct(...
        'modulation', 0.1, ...                      % 10% photon modulation against background
        'stimulusSamplingInterval',  1/10, ...      % 10 Hz stimulus refresh
        'responseTimeInterval', 5/1000, ...         % 5 milliseconds
        'photonNoise', true, ...
        'osNoise', false);
    end

    addCondition = true;
    if (addCondition)
    condData{numel(condData)+1} = struct(...
        'modulation', 0.1, ...                      % 10% photon modulation against background
        'stimulusSamplingInterval',  1/10, ...      % 10 Hz stimulus refresh
        'responseTimeInterval', 5/1000, ...         % 5 milliseconds
        'photonNoise', true, ...
        'osNoise', true);
    end

    addCondition = true;
    if (addCondition)
    condData{numel(condData)+1} = struct(...
        'modulation', 0.3, ...                      % 30% photon modulation against background
        'stimulusSamplingInterval',  1/10, ...      % 10 Hz stimulus refresh
        'responseTimeInterval', 5/1000, ...         % 5 milliseconds
        'photonNoise', true, ...
        'osNoise', true);
    end

end



% Run all the conditions
for stimulusConditionIndex = 1:numel(condData)
    % Get the condition data
    c = condData{stimulusConditionIndex};   
    
    % Run the simulation for this condition
    [theConeMosaic, theOIsequence, ...
        isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, ...
        stimulusTimeAxis, responseTimeAxis] = runSimulation(c.modulation, c.stimulusSamplingInterval, c.responseTimeInterval, c.photonNoise, c.osNoise);

    % Plot the results
    plotEverything(theConeMosaic, theOIsequence, isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, stimulusTimeAxis, responseTimeAxis, stimulusConditionIndex, c);
end

end

function [theConeMosaic, theOIsequence, ...
    isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, ...
    stimulusTimeAxis, responseTimeAxis] = runSimulation(modulation, stimulusSamplingInterval, responseTimeInterval, photonNoise, osNoise)

    % Define the time axis for the simulation (how much data we will generate)
    stimulusTimeAxis = -0.6:stimulusSamplingInterval:0.4;
    stimulusRampTau = 0.165;

    % Generate a uniform field scene with desired mean luminance
    FOV = 2.0; meanLuminance = 500;
    theScene = uniformFieldSceneCreate(FOV, meanLuminance);

    % Generate optics
    noOptics = false;
    theOI = oiGenerate(noOptics);

    % Generate the sequence of optical images
    theOIsequence = oiSequenceGenerateForRampedSceneModulation(theScene, theOI, stimulusTimeAxis, stimulusRampTau, modulation);

    % Generate the cone mosaic with eye movements for theOIsequence
    theConeMosaic = coneMosaicGenerate(photonNoise, osNoise, responseTimeInterval, stimulusSamplingInterval, numel(theOIsequence));

    % Compute mosaic response to sequence of OIs!
    [isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, responseTimeAxis] = ...
        computeMosaicResponse(theConeMosaic, theOIsequence, stimulusSamplingInterval, stimulusTimeAxis);

end


function  [isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, responseTimeAxis] = ...
    computeMosaicResponse(theConeMosaic, theOIsequence, stimulusSamplingInterval, stimulusTimeAxis)

    % Save a copy of the entire eye movement sequence
    eyeMovementsForOISequence = theConeMosaic.emPositions;

    % Check that this is possible
    eyeMovementFramesPerOpticalImage = stimulusSamplingInterval/theConeMosaic.os.timeStep;
    if (eyeMovementFramesPerOpticalImage < 1.0)
        error('\nEye movements per optical image: %f.\nEither decrease the responseTimeInterval or increase the stimulusSamplingInterval \n', eyeMovementFramesPerOpticalImage);
    end
    
    % Initialize our variables
    absorptions = []; lastEyeMovementIndex = 0; eyeMovementSequence = [];

    % Loop over the optical images and compute isomerizations
    for oiIndex = 1:numel(theOIsequence)

        % Retrieve eye movements for current OI
        firstEyeMovementIndex = lastEyeMovementIndex+1;
        lastEyeMovementIndex = round(oiIndex*eyeMovementFramesPerOpticalImage);
        eyeMovementIndices = firstEyeMovementIndex:lastEyeMovementIndex;
        eyeMovementPathForThisOI = eyeMovementsForOISequence(eyeMovementIndices,:);
        
        % Compute absorbions for current OI eye movement path
        theConeMosaic.emPositions = eyeMovementPathForThisOI;
        %fprintf('Added %d points for oiIndex: %d\n', lastEyeMovementIndex-firstEyeMovementIndex+1, oiIndex);
        absorptionsForThisOI = theConeMosaic.compute(theOIsequence{oiIndex}, 'currentFlag', false, 'newNoise', true);

        % Concatenate sequences
        absorptions = cat(3, absorptions, absorptionsForThisOI);
        eyeMovementSequence = cat(1, eyeMovementSequence, eyeMovementPathForThisOI);
    end

    % Define response time axis (starting at the origin of the stimulusTimeAxis)
    responseTimeAxis = stimulusTimeAxis(1) + (0:(size(absorptions,3)-1))*theConeMosaic.os.timeStep;

    % Convert to photon rate
    isomerizationRateSequence = absorptions/theConeMosaic.integrationTime;
    photoCurrentSequence = theConeMosaic.os.osCompute(isomerizationRateSequence, theConeMosaic.pattern, 'append', false);
end

function theConeMosaic = coneMosaicGenerate(photonNoise, osNoise, responseTimeInterval, stimulusSamplingInterval, opticalImageSequenceLength)
    % Generate a human cone mosaic with 1L, 1M and 1S cone
    theConeMosaic = coneMosaic;
    theConeMosaic.rows = 1;
    theConeMosaic.cols = 3;
    theConeMosaic.pattern = [2 3 4];
    theConeMosaic.noiseFlag = photonNoise;

    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = osNoise;
    % Set a custom timeStep, for @osLinear we do not need a very small value
    theOuterSegment.timeStep = responseTimeInterval;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;

    % Generate eye movement sequence for all oi's
    eyeMovementFramesPerOpticalImage = stimulusSamplingInterval/theConeMosaic.os.timeStep;
    theConeMosaic.emGenSequence(ceil(opticalImageSequenceLength*eyeMovementFramesPerOpticalImage));
end


function theOIsequence = oiSequenceGenerateForRampedSceneModulation(theScene, theOI, stimulusTimeAxis, stimulusRampTau, modulation)

    % Stimulus time ramp
    stimulusRamp = exp(-0.5*(stimulusTimeAxis/stimulusRampTau).^2);
    
    % Compute the optical image
    theOI = oiCompute(theOI, theScene);

    backgroundPhotons = oiGet(theOI, 'photons');
    for stimFrameIndex = 1:numel(stimulusTimeAxis)
        retinalPhotonsAtCurrentFrame = backgroundPhotons * (1.0 + modulation*stimulusRamp(stimFrameIndex));
        theOIsequence{stimFrameIndex} = oiSet(theOI, 'photons', retinalPhotonsAtCurrentFrame);
    end
end


function theOI = oiGenerate(noOptics)
    % Generate optics
    if (noOptics)
        theOI = oiCreate('diffraction limited');
        optics = oiGet(theOI,'optics');           
        optics = opticsSet(optics,'fnumber',0);
        optics = opticsSet(optics, 'off axis method', 'skip');
        theOI = oiSet(theOI,'optics', optics);
    else
        theOI = oiCreate('human');
    end
end

% Supporting functions
function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
    uniformScene = sceneCreate('uniformd65');
    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);
    % set the radiance (in photons/steradian/m^2/nm)
    photonFlux = 1e16;
    uniformScene = sceneSet(uniformScene, 'photons', photonFlux*ones(64,64,numel(sceneGet(uniformScene, 'wave'))));
    % set desired luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end

function plotEverything(theConeMosaic, theOIsequence, isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, stimulusTimeAxis, responseTimeAxis, figNo, condData)

    % Plot everything
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10+figNo*50 10+figNo*100 2300 560], 'Color', [1 1 1]);
    set(hFig, 'Name', sprintf('Modulation: %2.2f,     Stimulus Sampling: %2.1f ms,     Response Sampling: %2.1f ms,      PhotonNoise: %g,      osNoise: %g', condData.modulation, condData.stimulusSamplingInterval*1000, condData.responseTimeInterval*1000, condData.photonNoise, condData.osNoise));

    oiWavelengthAxis = oiGet(theOIsequence{1}, 'wave');
    
    %% Plot the photon rate at the center of the optical image
    subplot('Position', [0.03 0.07 0.18 0.89]);
    referencePositionOpticalImagePhotons = zeros(numel(oiWavelengthAxis), numel(theOIsequence));
    for oiIndex = 1:numel(theOIsequence)
        retinalPhotonsAtCurrentFrame = oiGet(theOIsequence{oiIndex}, 'photons');
        refRow = round(size(retinalPhotonsAtCurrentFrame,1)/2);
        refCol = round(size(retinalPhotonsAtCurrentFrame,2)/2);
        referencePositionOpticalImagePhotons(:, oiIndex) = squeeze(retinalPhotonsAtCurrentFrame(refRow, refCol, :));
    end
    hP = pcolor(stimulusTimeAxis, oiWavelengthAxis, referencePositionOpticalImagePhotons);
    set(hP, 'EdgeColor', 'none');
    hold on;
    % Plot the total photons (summed across all wavelengths)
    totalPhotons = sum(referencePositionOpticalImagePhotons,1);
    totalPhotonsNorm = oiWavelengthAxis(1) + (oiWavelengthAxis(end)-oiWavelengthAxis(1))*(totalPhotons-min(totalPhotons))/(max(totalPhotons)-min(totalPhotons));
    stairs(stimulusTimeAxis, totalPhotonsNorm, 'c-', 'LineWidth', 2.0);
    % Plot lines demarkating each OI time duration
    for oiIndex = 1:numel(theOIsequence)
        plot(stimulusTimeAxis(oiIndex)*[1 1], [oiWavelengthAxis(1) oiWavelengthAxis(end)], 'k-');
    end
    % Plot the origin in magenta
    plot([0 0], [oiWavelengthAxis(1) oiWavelengthAxis(end)], '-', 'Color', [0.7 0.1 0.3], 'LineWidth', 2);
    
    hold off; box on
    set(gca, 'YLim', [oiWavelengthAxis(1) oiWavelengthAxis(end)], 'XLim', [stimulusTimeAxis(1) stimulusTimeAxis(end)]);
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('wavelength (nm)', 'FontSize', 14, 'FontWeight', 'bold');
    title('optical image photon rate (image center)', 'FontSize', 14);
    hC = colorbar('Location', 'NorthOutside');
    hC.FontSize =  14;
    hC.Label.String = 'photons/sec';
    axis 'xy'
    colormap(gray(1024));

    %% Plot the eye movement sequence (different colors for different OIs)
    subplot('Position', [0.25 0.07 0.22 0.89]); hold on;
    eyeMovementRange = [-100 100];

    plot(responseTimeAxis, eyeMovementSequence(:,1)*theConeMosaic.pigment.width*1e6, '.', 'MarkerSize', 15, 'Color', 'r');
    hold on;
    plot(responseTimeAxis, eyeMovementSequence(:,2)*theConeMosaic.pigment.height*1e6, '.', 'MarkerSize', 15, 'Color', 'b');
    % Plot lines demarkating each OI time duration
    for oiIndex = 1:numel(theOIsequence)
        plot(stimulusTimeAxis(oiIndex)*[1 1], eyeMovementRange, 'k-');
    end
    % Plot the origin in magenta
    plot([0 0], eyeMovementRange, '-', 'Color', [0.7 0.1 0.3], 'LineWidth', 2);
    
    box on
    set(gca, 'YLim', [eyeMovementRange(1) eyeMovementRange(end)], 'XLim', [stimulusTimeAxis(1) stimulusTimeAxis(end)]);
    legend({'eye position (X)', 'eye position (Y)'});
    ylabel('X,Y eye position (microns)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('X,Y eye movements', 'FontSize', 14);
    
    
    %% Plot the LMS isomerizations
    subplot('Position', [0.50 0.07 0.22 0.89]);
    isomerizationRange = [min(isomerizationRateSequence(:)) 1.05*max(isomerizationRateSequence(:))];
    plot(responseTimeAxis, squeeze(isomerizationRateSequence(1,1,:)), 'r.', 'MarkerSize', 15, 'LineWidth', 1.5);
    hold on;
    plot(responseTimeAxis, squeeze(isomerizationRateSequence(1,2,:)), 'g.', 'MarkerSize', 15, 'LineWidth', 1.5);
    plot(responseTimeAxis, squeeze(isomerizationRateSequence(1,3,:)), 'b.', 'MarkerSize', 15, 'LineWidth', 1.5);
    
    % Plot lines demarkating each OI time duration
    for oiIndex = 1:numel(theOIsequence)
        plot(stimulusTimeAxis(oiIndex)*[1 1], [isomerizationRange(1) isomerizationRange(end)], 'k-');
    end
    % Plot the origin in magenta
    plot([0 0], isomerizationRange, '-', 'Color', [0.7 0.1 0.3], 'LineWidth', 2);
    
    hold off;
    set(gca, 'YLim', isomerizationRange, 'XLim', [stimulusTimeAxis(1) stimulusTimeAxis(end)]);
    ylabel('isomerization rate (R*/cone/sec)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('L,M,S-cone isomerization rates', 'FontSize', 14);

    % Plot the photocurrents
    subplot('Position', [0.75 0.07 0.22 0.89]);
    photoCurrentRange = [min(photoCurrentSequence(:)) max(photoCurrentSequence(:))+2];
    plot(responseTimeAxis, squeeze(photoCurrentSequence(1,1,:)), 'r.', 'MarkerSize', 15, 'LineWidth', 1.5);
    hold on;
    plot(responseTimeAxis, squeeze(photoCurrentSequence(1,2,:)), 'g.', 'MarkerSize', 15, 'LineWidth', 1.5);
    plot(responseTimeAxis, squeeze(photoCurrentSequence(1,3,:)), 'b.', 'MarkerSize', 15, 'LineWidth', 1.5);
    % Plot lines demarkating each OI time duration
    for oiIndex = 1:numel(theOIsequence)
        plot(stimulusTimeAxis(oiIndex)*[1 1], photoCurrentRange, 'k-');
    end
    % Plot the origin in magenta
    plot([0 0], photoCurrentRange, '-', 'Color', [0.7 0.1 0.3], 'LineWidth', 2);
    
    hold off;
    set(gca, 'XLim', [stimulusTimeAxis(1) stimulusTimeAxis(end)], 'YLim', photoCurrentRange);
    ylabel('photocurrent (pA)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('@osLinear response', 'FontSize', 14);
    
    drawnow
end

