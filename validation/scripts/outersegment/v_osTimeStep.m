function varargout = v_osTimeStep(varargin)
%
% Demonstrate simulations using three different timebases, one for stimuli (based on stimulus refresh rate), 
% one for absorptions and eye movements (based on coneMosaic.integrationTime), and a third one for 
% outer segment current computations (based on os.timeStep) 
% Also demonstrates computations based on multi-optical images with eye mvoements.
%
% NPC, ISETBIO TEAM, 2016
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;

% Examine effects of varying the response time interval (os.timeStep)
conditionSet = 1;

% Examine the effects of varying the integrationTime
%conditionSet = 2;

% Examine the magnitudes of the photon and OS noise
conditionSet = 3;

condData = makeConditionSet(conditionSet);

% Run all the conditions
for stimulusConditionIndex = 1:numel(condData)
    % Get the condition data
    c = condData{stimulusConditionIndex};   
    
    % Run the simulation for this condition
    tic
    [theConeMosaic, theOIsequence, ...
        isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, ...
        stimulusTimeAxis, absorptionsTimeAxis, responseTimeAxis] = runSimulation(c.mosaicSize, c.meanLuminance, c.modulation, c.modulationRegion, c.stimulusSamplingInterval, c.integrationTime, c.responseTimeInterval, c.photonNoise, c.osNoise);
    toc
    
    % Plot the results
    plotEverything(theConeMosaic, theOIsequence, isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, stimulusTimeAxis, absorptionsTimeAxis, responseTimeAxis, stimulusConditionIndex, c);
end

end

function [theConeMosaic, theOIsequence, ...
    isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, ...
    stimulusTimeAxis, absorptionsTimeAxis, responseTimeAxis] = runSimulation(mosaicSize, meanLuminance, modulation, modulationRegion, stimulusSamplingInterval, integrationTime, responseTimeInterval, photonNoise, osNoise)

    % Define the time axis for the simulation (how much data we will generate)
    stimulusTimeAxis = -0.6:stimulusSamplingInterval:0.4;
    stimulusRampTau = 0.165;

    % Generate a uniform field scene with desired mean luminance
    if (isnan(mosaicSize))
        FOV = 0.5;
    else
        FOV = max(mosaicSize);
    end
    theScene = uniformFieldSceneCreate(FOV, meanLuminance);

    % Generate optics
    noOptics = false;
    theOI = oiGenerate(noOptics);

    % Generate the sequence of optical images
    theOIsequence = oiSequenceGenerateForRampedSceneModulation(theScene, theOI, stimulusTimeAxis, stimulusRampTau, modulation, modulationRegion);

    % Generate the cone mosaic with eye movements for theOIsequence
    theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, responseTimeInterval, stimulusSamplingInterval, numel(theOIsequence));

    % Compute mosaic response to sequence of OIs!
    [isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, absorptionsTimeAxis, responseTimeAxis] = ...
        computeMultiOIMosaicResponse(theConeMosaic, theOIsequence, stimulusSamplingInterval, stimulusTimeAxis);
end


function  [isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, absorptionsTimeAxis, responseTimeAxis] = ...
    computeMultiOIMosaicResponse(theConeMosaic, theOIsequence, stimulusSamplingInterval, stimulusTimeAxis)

    % Save a copy of the entire eye movement sequence
    eyeMovementsForOISequence = theConeMosaic.emPositions;

    % Check that there is at least 1 eye movement per OI
    eyeMovementFramesPerOpticalImage = stimulusSamplingInterval/theConeMosaic.integrationTime;
    if (eyeMovementFramesPerOpticalImage < 1.0)
        error('\nEye movements per optical image: %f.\nEither decrease the responseTimeInterval or increase the stimulusSamplingInterval \n', eyeMovementFramesPerOpticalImage);
    end
    
    % Initialize our time series
    absorptions = []; absorptionsTimeAxis = []; eyeMovementSequence = [];  
    lastEyeMovementIndex = 0;
    
    % Loop over the optical images and compute isomerizations
    for oiIndex = 1:numel(theOIsequence)

        % Retrieve eye movements for current OI
        firstEyeMovementIndex = lastEyeMovementIndex+1;
        lastEyeMovementIndex = round(oiIndex*eyeMovementFramesPerOpticalImage);
        eyeMovementIndices = firstEyeMovementIndex:lastEyeMovementIndex;
        eyeMovementPathForThisOI = eyeMovementsForOISequence(eyeMovementIndices,:);
        
        % Compute absorptions for current OI and eye movement path
        absorptionsForThisOI = theConeMosaic.compute(theOIsequence{oiIndex}, ...
            'emPath', eyeMovementPathForThisOI, ...      % current OI eye movement path
            'currentFlag', false, ...
            'newNoise', true);

        % Concatenate sequences
        if (isempty(absorptionsTimeAxis))
            absorptionsTimeAxis = theConeMosaic.absorptionsTimeAxis;
            eyeMovementSequence = eyeMovementPathForThisOI;
            absorptions = absorptionsForThisOI;
        else
            absorptionsTimeAxis = cat(2, absorptionsTimeAxis, absorptionsTimeAxis(end)+theConeMosaic.absorptionsTimeAxis);
            absorptions = cat(3, absorptions, absorptionsForThisOI);
            eyeMovementSequence = cat(1, eyeMovementSequence, eyeMovementPathForThisOI);
        end
    end

    % isomerization rate sequence in integrationTime time units
    
    isomerizationRateSequence = absorptions / theConeMosaic.integrationTime;
    
    % compute the os time axis
    dtOS = theConeMosaic.os.timeStep;
    osTimeAxis = absorptionsTimeAxis(1): dtOS :absorptionsTimeAxis(end);
    
    % Resample absorptions to osTimeAxis
    resampledAbsorptionsSequence = coneMosaic.resampleAbsorptionsSequence(absorptions, absorptionsTimeAxis, osTimeAxis);
    
    % Convert to photon rate
    pRate = resampledAbsorptionsSequence/dtOS;
    photoCurrentSequence = theConeMosaic.os.osCompute(pRate, theConeMosaic.pattern, 'append', false);
    
    % Define response time axis (starting at the origin of the stimulusTimeAxis)
    responseTimeAxis = stimulusTimeAxis(1) + osTimeAxis;
    absorptionsTimeAxis = stimulusTimeAxis(1) + absorptionsTimeAxis;    
end

function theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, responseTimeInterval, stimulusSamplingInterval, opticalImageSequenceLength)
    
    % Default human mosaic
    theConeMosaic = coneMosaic;
    
    % Adjust size
    if isnan(mosaicSize)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicSize);
    end
    
    % Set the noise
    theConeMosaic.noiseFlag = photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = integrationTime;
    
    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = osNoise;
    % Set a custom timeStep, for @osLinear we do not need a very small value
    theOuterSegment.timeStep = responseTimeInterval;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;

    % Generate eye movement sequence for all oi's
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/theConeMosaic.integrationTime;
    if (eyeMovementsNumPerOpticalImage == floor(eyeMovementsNumPerOpticalImage))
        fprintf('Will generate %d eye movements\n', eyeMovementsNumPerOpticalImage*opticalImageSequenceLength)
        theConeMosaic.emGenSequence(eyeMovementsNumPerOpticalImage*opticalImageSequenceLength);
    else
       error(' eyeMovementsNumPerOpticalImage (%g) is not integer.\nStimulus sampling interval:%g integration time: %g\n', eyeMovementsNumPerOpticalImage, stimulusSamplingInterval, theConeMosaic.integrationTime);
    end
    % * * * * * * SUPER IMPORTANT * * * * *
    % Set the integration time correctly, because emGenSequence re-sets it to os.timeStep * eyeMovementFrames
    theConeMosaic.integrationTime
end


function theOIsequence = oiSequenceGenerateForRampedSceneModulation(theScene, theOI, stimulusTimeAxis, stimulusRampTau, modulation, modulationRegion)
    % Stimulus time ramp
    stimulusRamp = exp(-0.5*(stimulusTimeAxis/stimulusRampTau).^2);
    
    % Compute the optical image
    theOI = oiCompute(theOI, theScene);
    backgroundPhotons = oiGet(theOI, 'photons');
    
    fprintf('Computing sequence of optical images\n');
    
    for stimFrameIndex = 1:numel(stimulusTimeAxis)
        if strcmp(modulationRegion, 'FULL')
            retinalPhotonsAtCurrentFrame = backgroundPhotons * (1.0 + modulation*stimulusRamp(stimFrameIndex));
        elseif strcmp(modulationRegion, 'CENTER')
            if (stimFrameIndex == 1)
                pos = oiGet(theOI, 'spatial support', 'microns');
                ecc = sqrt(squeeze(sum(pos.^2, 3)));
                idx = find(ecc < 0.25*max(pos(:)));
                [idx1, idx2] = ind2sub(size(ecc), idx);
            end
            retinalPhotonsAtCurrentFrame = backgroundPhotons;
            retinalPhotonsAtCurrentFrame(idx1, idx2, :) = retinalPhotonsAtCurrentFrame(idx1, idx2, :) * (1.0 + modulation*stimulusRamp(stimFrameIndex));
        else
            error('Unknown modulationRegion ''%s'', modulationRegion');
        end
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


function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
    uniformScene = sceneCreate('uniformd65');
    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);
    % set some radiance (in photons/steradian/m^2/nm)
    photonFlux = 1e16;
    uniformScene = sceneSet(uniformScene, 'photons', photonFlux*ones(64,64,numel(sceneGet(uniformScene, 'wave'))));
    % adjust radiance according to desired  mean luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end


function condData = makeConditionSet(conditionSet)

    % scene mean luminance
    meanLuminance = 1500;    
    
    % Assemble conditions to run.
    condData = {};

    switch conditionSet

        % Examine effects of varying the responseTimeInterval (os.timeStep)
        case 1

            % Steady params
            c0 = struct(...
                'mosaicSize', nan, ...                      % 1 L-, 1 M-, and 1 S-cone only
                'meanLuminance', meanLuminance, ...         % scene mean luminance
                'modulation', 0.5, ...                      % 50%  modulation against background
                'modulationRegion', 'FULL', ...             % modulate the entire image (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/10, ...      % 100 Hz stimulus refresh
                'integrationTime', 20/1000, ...             % 20 milliseconds
                'responseTimeInterval', nan, ...            % we'll vary that one
                'photonNoise', true, ...
                'osNoise', false);

            % Varied params
            c0.responseTimeInterval = 20/1000;              % 20 millisecond response interval
            condData{numel(condData)+1} = c0;
            
            c0.responseTimeInterval = 10/1000;              % 10 millisecond response interval
            condData{numel(condData)+1} = c0;
            
            c0.responseTimeInterval = 2/1000;               % 5 millisecond response interval
            condData{numel(condData)+1} = c0;
            

        % Effects of varying the integrationTime
        case 2

            % Steady params
            c0 = struct(...
                'mosaicSize', nan, ...                      % 1 L-, 1 M-, and 1 S-cone only
                'meanLuminance', meanLuminance, ...         % scene mean luminance
                'modulation', 0.5, ...                      % 50%  modulation against background
                'modulationRegion', 'FULL', ...             % modulate the entire image (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/10, ...      % 10 Hz stimulus refresh
                'responseTimeInterval', 1/1000, ...         % 1 milliseconds
                'integrationTime', nan, ...                 % we will vary this one
                'photonNoise', true, ...
                'osNoise', false);
            
            % Varied params
            c0.integrationTime = 100/1000;                          % 100 ms
            condData{numel(condData)+1} = c0;
            
            c0.integrationTime = 10/1000;                           % 10 ms
            condData{numel(condData)+1} = c0;
            
            c0.integrationTime = 1/1000;                            % 1 ms
            condData{numel(condData)+1} = c0;
            
            

        case 3    
            
            stimulusRefreshRateInHz = 25;
            eyeMovementsPerStimulusRefresh = 6;
            
            % Steady params
            c0 = struct(...
                'mosaicSize', nan, ...                      % 1 L-, 1 M-, and 1 S-cone only
                'meanLuminance', meanLuminance, ...         % scene mean luminance
                'modulation', 0.5, ...                      % % modulation against background
                'modulationRegion', 'CENTER', ...           % modulate the center only  (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/stimulusRefreshRateInHz, ...      % 87 Hz stimulus refresh
                'integrationTime', 5/1000, ...             % 5 milliseconds
                'responseTimeInterval', 1/(stimulusRefreshRateInHz*eyeMovementsPerStimulusRefresh), ...  
                'photonNoise', nan, ...
                'osNoise', nan);
            
            % Varied params
            % No noise
            c0.photonNoise = false;
            c0.osNoise = false;
            condData{numel(condData)+1} = c0;
            
            % Photon noise
            c0.photonNoise = true;
            c0.osNoise = false;
            condData{numel(condData)+1} = c0;
            
            % OS noise 
            c0.photonNoise = false;
            c0.osNoise = true;
            condData{numel(condData)+1} = c0;
            
            % Both photon noise and OS noise
            c0.photonNoise = true;
            c0.osNoise = true;
            condData{numel(condData)+1} = c0;
            
        case 4
            
            % Custom
            stimulusRefreshRateInHz = 30;
            eyeMovementsPerStimulusRefresh = 6;
            
            % Steady params
            c0 = struct(...
                'mosaicSize', 0.2, ...               % FOV = 3x3 deg
                'meanLuminance', meanLuminance, ...         % scene mean luminance
                'modulation', 0.5, ...                      % % modulation against background
                'modulationRegion', 'CENTER', ...           % modulate the center only  (choose b/n 'FULL', and 'CENTER')
                'stimulusSamplingInterval',  1/stimulusRefreshRateInHz, ...      % 87 Hz stimulus refresh
                'integrationTime', 20/1000, ...             % 20 milliseconds
                'responseTimeInterval', 1/(stimulusRefreshRateInHz*eyeMovementsPerStimulusRefresh), ...       % 2 eye movements / stimulus frame
                'photonNoise', nan, ...
                'osNoise', nan);
            
    end
end

function plotEverything(theConeMosaic, theOIsequence, isomerizationRateSequence, photoCurrentSequence, eyeMovementSequence, stimulusTimeAxis, absorptionsTimeAxis, responseTimeAxis, figNo, condData)

    % Plot the sequence of OIs with the eye movements 
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10+figNo*50 10+figNo*100 1600 800], 'Color', [1 1 1]);
    set(hFig, 'Name', sprintf('Scene Mean Luminance: %2.1f cd/m2,     Modulation: %2.2f,     Stimulus Sampling: %2.1f ms,    integrationTime: %2.1f ms,    OS Response Sampling: %2.1f ms', condData.meanLuminance, condData.modulation, condData.stimulusSamplingInterval*1000, condData.integrationTime*1000, condData.responseTimeInterval*1000));

    plotRows = round(0.75*sqrt(numel(theOIsequence))); 
    plotCols = ceil(numel(theOIsequence)/plotRows);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', plotRows, ...
           'colsNum', plotCols, ...
           'heightMargin',   0.03, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02);
       
    maxRGB = 0;
    for oiIndex = 1:numel(theOIsequence)
        tmp = xyz2rgb(oiGet(theOIsequence{oiIndex}, 'xyz'));
        if (maxRGB < max(tmp(:)))
            maxRGB = max(tmp(:));
        end
         oiImage{oiIndex} = tmp;
    end
    
    for oiIndex = 1:numel(theOIsequence)
        pos = oiGet(theOIsequence{oiIndex}, 'spatial support', 'microns');
        oiXaxis = pos(1,:,1); oiYaxis = pos(:,1,2);
        r = floor((oiIndex-1)/plotCols)+1;
        c = mod((oiIndex-1), plotCols)+1;
        subplot('Position', subplotPosVectors(r,c).v);
        
        % Plot the OI
        imagesc(oiXaxis, oiYaxis, oiImage{oiIndex}/maxRGB);
        % Overlay the eye movements up to this point
        hold on;
        idx = find(absorptionsTimeAxis <= stimulusTimeAxis(oiIndex));
        plot(eyeMovementSequence(idx,1)*theConeMosaic.pigment.width*1e6, eyeMovementSequence(idx,2)*theConeMosaic.pigment.width*1e6, 'k.-', 'LineWidth', 1.0);
        if (oiIndex < numel(theOIsequence))
            idx = find((absorptionsTimeAxis>=stimulusTimeAxis(oiIndex)) & (absorptionsTimeAxis<stimulusTimeAxis(oiIndex+1)));
        else
            idx = find((absorptionsTimeAxis>=stimulusTimeAxis(oiIndex)));
        end
        % Emphasize in red, the eye movements for the current framer
        plot(eyeMovementSequence(idx,1)*theConeMosaic.pigment.width*1e6, eyeMovementSequence(idx,2)*theConeMosaic.pigment.width*1e6, 'r.-', 'LineWidth', 4.0);
        
        % overlay the cone mosaic
        if (oiIndex == 1)
            plot(theConeMosaic.coneLocs(:,1)*1e6, theConeMosaic.coneLocs(:,2)*1e6, 'k.');
        end
        
        axis 'image'; axis 'xy';
        set(gca, 'CLim', [0 1]);
        if ~((r == plotRows) && (c == 1))
            set(gca, 'XTick', [], 'YTick', []);
        end
        title(sprintf('t: %2.2f', stimulusTimeAxis(oiIndex)));
    end
    
    % Plot time-varying responses
    hFig = figure(figNo+1000); clf;
    set(hFig, 'Position', [10+figNo*50 10+figNo*100 2300 560], 'Color', [1 1 1]);
    set(hFig, 'Name', sprintf('Scene Mean Luminance: %2.1f cd/m2,     Modulation: %2.2f,     Stimulus Sampling: %2.1f ms,     Integration Time: %2.1f ms,   Response Sampling: %2.1f ms,      PhotonNoise: %g,      osNoise: %g', condData.meanLuminance, condData.modulation, condData.stimulusSamplingInterval*1000, condData.integrationTime*1000, condData.responseTimeInterval*1000, condData.photonNoise, condData.osNoise));

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

    plot(absorptionsTimeAxis, eyeMovementSequence(:,1)*theConeMosaic.pigment.width*1e6, '.', 'MarkerSize', 15, 'Color', 'r');
    hold on;
    plot(absorptionsTimeAxis, eyeMovementSequence(:,2)*theConeMosaic.pigment.height*1e6, '.', 'MarkerSize', 15, 'Color', 'b');
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
    if (theConeMosaic.rows ==1) && (theConeMosaic.cols == 3)
       referenceConeRows = [1 1 1]; referenceConeCols = [1 2 3];
    else
       % Find the (row,col) coords of the center-most L, M and S-cone
       for k = 1:3
            coneIndices = find(theConeMosaic.pattern == k+1); 
            [~, idx] = min(sum((theConeMosaic.coneLocs(coneIndices,:)).^2, 2));
            [referenceConeRows(k), referenceConeCols(k)] = ind2sub(size(theConeMosaic.pattern), coneIndices(idx));
       end
    end
    subplot('Position', [0.50 0.07 0.22 0.89]);
    isomerizationRange = [min(isomerizationRateSequence(:)) 1.05*max(isomerizationRateSequence(:))];
    hold  on
    coneColors = [1 0 0; 0 1 0; 0 0 1];
    for k = 1:3
        plot(absorptionsTimeAxis, squeeze(isomerizationRateSequence(referenceConeRows(k),referenceConeCols(k),:)), '.', 'Color', squeeze(coneColors(k,:)), 'MarkerSize', 15, 'LineWidth', 1.5);
    end
    
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

    %% Plot the photocurrents
    subplot('Position', [0.75 0.07 0.22 0.89]);
    photoCurrentRange = [min(photoCurrentSequence(:)) max(photoCurrentSequence(:))+2];
    hold on;
    for k = 1:3
        plot(responseTimeAxis, squeeze(photoCurrentSequence(referenceConeRows(k),referenceConeCols(k),:)), 'k.', 'Color', squeeze(coneColors(k,:)), 'MarkerSize', 15, 'LineWidth', 1.5);
    end
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

