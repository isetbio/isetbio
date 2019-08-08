function t_osTimeStep
% Show how to change timebase of os calculation.
%
% Syntax:
%   t_osTimeStep
%
% Description:
%    Demonstrate simulations using three different timebases, one for
%    stimuli (based on stimulus refresh rate), one for absorptions and eye
%    movements (based on coneMosaic.integrationTime), and a third one for
%    outer segment current computations (based on os.timeStep)
%
%    Also demonstrates usage of the computeForOISequence() method of
%    @coneMosaic, which computes absorptions and photocurrents for a
%    sequence of sequentially presented optical images with eye movements.
%
%    This is a rather complex tutorial that can demonstrate lots of things.
%    It has many subfunctions which do lots of the work, and which have no
%    comments.
%
%    [Note: DHB - This needs more extensive comments before it can be
%    used in an introductory way.]
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    XX/XX/16  NPC  ISETBIO TEAM, 2016
%    08/07/17  dhb  Commenting pass. Moved to Advanced.
%                   Decided this was too hard and moved to advanced.
%    11/19/18  JNM  Formatting

%% Init
ieInit;

%% Examine the effects of varying the integrationTime
conditionSet = 1;

% % Examine effects of varying the response time interval (os.timeStep)
% conditionSet = 2;

% % Examine and contrast the magnitudes of photon noise vs OS noise
% conditionSet = 3;

% % Custom condition
% conditionSet = 4;

condData = makeConditionSet(conditionSet);

% Run all the conditions
for stimulusConditionIndex = 1:numel(condData)
    % Get the condition data
    c = condData{stimulusConditionIndex};

    % Run the simulation for this condition
    [theConeMosaic, theOIsequence, ...
        isomerizationRateSequence, photoCurrentSequence, ...
        oiTimeAxis, absorptionsTimeAxis, photoCurrentTimeAxis] = ...
        runSimulation(c.mosaicSize, c.meanLuminance, c.modulation, ...
        c.modulationRegion, c.stimulusSamplingInterval, ...
        c.integrationTime, c.osTimeStep, c.photonNoise, c.osNoise);

    % Plot the results
    plotEverything(theConeMosaic, theOIsequence, ...
        isomerizationRateSequence, photoCurrentSequence, oiTimeAxis, ...
        absorptionsTimeAxis, photoCurrentTimeAxis, ...
        stimulusConditionIndex, c);
end

end

function [theConeMosaic, theOIsequence, ...
    isomerizationRateSequence, photoCurrentSequence, ...
    oiTimeAxis, absorptionsTimeAxis, photoCurrentTimeAxis] = ...
    runSimulation(mosaicSize, meanLuminance, modulationGain, ...
    modulationRegion, stimulusSamplingInterval, integrationTime, ...
    osTimeStep, photonNoise, osNoise)
% Embedded helper function to run the simulation
%
% Syntax:
%   [theConeMosaic, theOIsequence, ...
%       isomerizationRateSequence, photoCurrentSequence, ...
%       oiTimeAxis, absorptionsTimeAxis, photoCurrentTimeAxis] = ...
%       runSimulation(mosaicSize, meanLuminance, modulationGain, ...
%       modulationRegion, stimulusSamplingInterval, integrationTime, ...
%       osTimeStep, photonNoise, osNoise)
%
% Description:
%    An embedded helper function in t_osTimeStep designed to run the
%    desired simulation.
%
% Inputs:
%    mosaicSize                - Vector. A vector of the mosaic size. Can
%                                also contain 'NaN' if want to use default.
%                                Default assumes FOV of 0.2.
%    meanLuminance             - Numeric. The mean luminance.
%    modulationGain            - Numeric. The modulation gain.
%    modulationRegion          - String. The string describing the
%                                modulation region.
%    stimulusSamplingInterval  - Numeric. The stimulus sampling interval.
%    integrationTime           - Numeric. The integration time.
%    osTimeStep                - Numeric. The outersegment time step.
%    photonNoise               - String. A string describing the type of
%                                photon noise.
%    osNoise                   - String. A string describing the
%                                outersegment noise type.
%
% Outputs:
%    theConeMosaic             - Object. A created cone mosaic object.
%    theOIsequence             - Object. A created OI sequence object.
%    isomerizationRateSequence - Matrix. A 4-D numeric matrix containing
%                                the isomerization rate sequence.
%    photoCurrentSequence      - Matrix. A 4-D numeric matrix containing
%                                the photocurrent sequence.
%    oiTimeAxis                - Vector. A vector of the OI time axis
%    absorptionsTimeAxis       - Vector. A vector containing the
%                                absorptions time axis.
%    photoCurrentTimeAxis      - Vector. A vector containing the
%                                photocurrent time axis.
%
% Optional key/value pairs:
%    None.
%

    % Define the time axis for the simulation
    minTime = -0.84;
    maxTime = 0.72;
    oiTimeAxis = minTime:stimulusSamplingInterval:maxTime;

    % Compute the stimulus modulation function
    stimulusRampTau = 0.18;
    modulationFunction = modulationGain * ...
        exp(-0.5 * (oiTimeAxis / stimulusRampTau) .^ 2);

    % Generate a uniform field scene with desired mean luminance
    if (isnan(mosaicSize))
        FOV = 0.2;
    else
        FOV = max(mosaicSize);
    end
    theScene = uniformFieldSceneCreate(FOV, meanLuminance);

    % Generate optics
    noOptics = false;
    theOI = oiGenerate(noOptics);

    % Generate the sequence of optical images
    theOIsequence = oiSequenceGenerate(theScene, theOI, oiTimeAxis, ...
        modulationFunction, modulationRegion);

    % Generate the cone mosaic with eye movements for theOIsequence
    theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, ...
        osNoise, integrationTime, osTimeStep);

    % Generate eye movement sequence for all oi's
    eyeMovementsNum = ...
        theOIsequence.maxEyeMovementsNumGivenIntegrationTime(...
        theConeMosaic.integrationTime);
    theConeMosaic.emGenSequence(eyeMovementsNum);

    [absorptionsCountSequence, photoCurrentSequence] = ...
            theConeMosaic.computeForOISequence(theOIsequence, ...
            'currentFlag', true);
    absorptionsTimeAxis = theConeMosaic.timeAxis ...
        + theOIsequence.timeAxis(1);
    photoCurrentTimeAxis = absorptionsTimeAxis;

    %theOIsequence.visualizeWithEyeMovementSequence(absorptionsTimeAxis);

    % Compute photon rate from photon count
    isomerizationRateSequence = absorptionsCountSequence ...
        / theConeMosaic.integrationTime;
end

% ------- Helper functions --------
function theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, ...
    osNoise, integrationTime, osTimeStep)
% Embedded function to gerneate a human cone mosaic.
%
% Syntax:
%   theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, ...
%       osNoise, integrationTime, osTimeStep)
%
% Description:
%    An embedded helper function inside of t_osTimeStep, designed to create
%    a human cone mosaic.
%
% Inputs:
%    mosaicSize      - Vector. A vector containing the mosaic size. To use
%                      default, enter NaN.
%    photonNoise     - String. A string describing the photon noise type.
%    osNoise         - String. A string describing the type of OS noise.
%    integrationTime - Numeric. The integration time.
%    osTimeStep      - Numeric. The OS time step.
%
% Outputs:
%    theConeMosaic   - Object. A created cone mosaic object.
%
% Optional key/value pairs:
%    None.
%

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

    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = osTimeStep;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;
end

function theOIsequence = oiSequenceGenerate(theScene, theOI, ...
    oiTimeAxis, modulationFunction, modulationType)
% Helper function to generate OI Sequence
%
% Syntax:
%   theOISequence = oiSequenceGenerate(theScene, theOI, oiTimeAxis, ...
%       modulationFunction, modulationType)
%
% Description:
%    An embedded helper function in t_osTimeStep designed to generate the
%    necessary OI Sequence for the function.
%
% Inputs:
%    theScene           - Struct. A scene structure.
%    theOI              - Struct. An optical image structure.
%    oiTimeAxis         - Vector. A vector describing the OI time axis.
%    modulationFunction - Vector. A vector containing the modulation
%                         function for the OI sequence.
%    modulationType     - String. A string describing the modulation type.
%
% Outputs:
%    theOISequence      - Object. An optical image sequence object.
%
% Optional key/value pairs:
%    None.
%

    % Compute the background and modulated optical images
    oiBackground = oiCompute(theOI, theScene);
    oiModulated  = oiBackground;

    if strcmp(modulationType, 'FULL')
        theOIsequence = oiSequence(oiBackground, oiModulated, ...
            oiTimeAxis, modulationFunction, 'composition', 'add');
    else
        pos = oiGet(oiBackground, 'spatial support', 'microns');
        modulationRegion.radiusInMicrons = 0.75 * max(pos(:));
        theOIsequence = oiSequence(oiBackground, oiModulated, ...
            oiTimeAxis, modulationFunction, ...
            'modulationRegion', modulationRegion);
    end
end

function theOI = oiGenerate(noOptics)
    % Generate optics
    if (noOptics)
        theOI = oiCreate('diffraction limited');
        optics = oiGet(theOI, 'optics');
        optics = opticsSet(optics, 'fnumber', 0);
        optics = opticsSet(optics, 'off axis method', 'skip');
        theOI = oiSet(theOI, 'optics', optics);
    else
        theOI = oiCreate('human');
    end
end

function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
% A helper function to create a scene with an uniform field
%
% Syntax:
%   uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
%
% Description:
%    An embedded helper function in t_osTimeStep designed to create a scene
%    containing an uniform field.
%
% Inputs:
%    FOV           - Numeric. The field of view.
%    meanLuminance - Numeric. The mean luminance.
%
% Outputs:
%    uniformScene  - Struct. A created scene structure.
%
% Optional key/value pairs:
%    None.
%
    uniformScene = sceneCreate('uniform equal photon', 128);
    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);
    % adjust radiance according to desired  mean luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end

function plotEverything(theConeMosaic, theOIsequence, ...
    isomerizationRateSequence, photoCurrentSequence, oiTimeAxis, ...
    absorptionsTimeAxis, responseTimeAxis, figNo, condData)

    % Plot the sequence of OIs with the eye movements
    hFig = figure(figNo); clf;
    set(hFig, ...
        'Position', [10 + figNo * 50, 10 + figNo * 100, 1920, 760], ...
        'Color', [1 1 1]);
    set(hFig, 'Name', sprintf(strcat("Scene Mean Luminance: %2.1f ", ...
        "cd/m2, Modulation: %2.2f, Stimulus Sampling: %2.1f ", ...
        "ms, Integration Time: %2.1f ms, osTimeStep: %2.1f ms, ", ...
        "      PhotonNoise: %s, osNoise: %s"), ...
        condData.meanLuminance, condData.modulation, ...
        condData.stimulusSamplingInterval * 1000, ...
        condData.integrationTime * 1000, condData.osTimeStep * 1000, ...
        condData.photonNoise, condData.osNoise));

    tabGroup = uitabgroup('Parent', hFig);

    eyeMovementsTab = uitab(tabGroup, ...
        'Title', '2D eye movements and OI sequence', ...
        'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 1]);
    timeSeriesTab = uitab(tabGroup, ...
        'Title', strcat("Time series: optical image photon rate, ", ...
        "eye movements, aborptions, & photocurrents"), ...
        'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 1]);

    set(tabGroup, 'SelectedTab', eyeMovementsTab);
    axes('parent', eyeMovementsTab);

    plotRows = round(0.75 * sqrt(theOIsequence.length));
    plotCols = ceil(theOIsequence.length/plotRows);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', plotRows, ...
           'colsNum', plotCols, ...
           'heightMargin', 0.03, ...
           'widthMargin', 0.03, ...
           'leftMargin', 0.03, ...
           'rightMargin', 0.00, ...
           'bottomMargin', 0.03, ...
           'topMargin', 0.02);

    maxRGB = 0;
    for oiIndex = 1:theOIsequence.length
        tmp = xyz2rgb(oiGet(theOIsequence.frameAtIndex(oiIndex), 'xyz'));
        if (maxRGB < max(tmp(:)))
            maxRGB = max(tmp(:));
        end
         oiImage{oiIndex} = tmp;
    end

    % Retrieve the eye movement sequence
    eyeMovementSequence = theConeMosaic.emPositions;

    for oiIndex = 1:theOIsequence.length
        pos = oiGet(theOIsequence.frameAtIndex(oiIndex), ...
            'spatial support', 'microns');
        oiXaxis = pos(1, :, 1);
        oiYaxis = pos(:, 1, 2);
        r = floor((oiIndex - 1) / plotCols) + 1;
        c = mod((oiIndex - 1), plotCols) + 1;
        subplot('Position', subplotPosVectors(r, c).v);

        % Plot the OI
        imagesc(oiXaxis, oiYaxis, oiImage{oiIndex}/maxRGB);
        % Overlay the eye movements up to this point
        hold on;
        idx = find(absorptionsTimeAxis < oiTimeAxis(oiIndex));
        plot(eyeMovementSequence(idx, 1) * theConeMosaic.pigment.width ...
            * 1e6, eyeMovementSequence(idx, 2) ...
            * theConeMosaic.pigment.width * 1e6, 'ks-', ...
            'LineWidth', 1.5, 'MarkerFaceColor', [0.3 0.3 0.3]);
        if (oiIndex < theOIsequence.length)
            idx = find((absorptionsTimeAxis >= oiTimeAxis(oiIndex)) ...
                & (absorptionsTimeAxis<oiTimeAxis(oiIndex+1)));
        else
            idx = find((absorptionsTimeAxis >= oiTimeAxis(oiIndex)));
        end
        % Emphasize in red, the eye movements for the current framer
        plot(eyeMovementSequence(idx, 1) * theConeMosaic.pigment.width ...
            * 1e6, eyeMovementSequence(idx, 2) ...
            * theConeMosaic.pigment.width * 1e6, 'rs-', ...
            'LineWidth', 1.5, 'MarkerFaceColor', [1.0 0.5 0.5]);

        % overlay the cone mosaic
        if (oiIndex == 1)
            plot(theConeMosaic.coneLocs(:, 1) * 1e6, ...
                theConeMosaic.coneLocs(:, 2) * 1e6, 'k.');
        end

        axis 'image';
        axis 'xy';
        set(gca, 'CLim', [0 1], 'FontSize', 14);
        if ~((r == plotRows) && (c == 1))
            set(gca, 'XTick', [], 'YTick', []);
        end
        title(sprintf('t: %2.1f msec', 1000 * oiTimeAxis(oiIndex)), ...
            'FontSize', 14);
    end

    % Plot time-series responses
    set(tabGroup, 'SelectedTab', timeSeriesTab);
    axes('parent', timeSeriesTab);

    %% Plot the photon rate at the center of the optical image
    subplot('Position', [0.03 0.07 0.18 0.89]);
    oiWavelengthAxis = oiGet(theOIsequence.frameAtIndex(1), 'wave');
    referencePositionOpticalImagePhotons = ...
        zeros(numel(oiWavelengthAxis), theOIsequence.length);
    for oiIndex = 1:theOIsequence.length
        retinalPhotonsAtCurrentFrame = ...
            oiGet(theOIsequence.frameAtIndex(oiIndex), 'photons');
        refRow = round(size(retinalPhotonsAtCurrentFrame, 1) / 2);
        refCol = round(size(retinalPhotonsAtCurrentFrame, 2) / 2);
        referencePositionOpticalImagePhotons(:, oiIndex) = ...
            squeeze(retinalPhotonsAtCurrentFrame(refRow, refCol, :));
    end
    hP = pcolor(oiTimeAxis, oiWavelengthAxis, ...
        referencePositionOpticalImagePhotons);
    set(hP, 'EdgeColor', 'none');
    hold on;
    % Plot the total photons (summed across all wavelengths)
    totalPhotons = sum(referencePositionOpticalImagePhotons, 1);
    totalPhotonsNorm = oiWavelengthAxis(1) + (oiWavelengthAxis(end) ...
        - oiWavelengthAxis(1)) * (totalPhotons - min(totalPhotons)) ...
        / (max(totalPhotons) - min(totalPhotons));
    stairs(oiTimeAxis, totalPhotonsNorm, 'c-', 'LineWidth', 2.0);
    plotStimulusTimes([min(totalPhotonsNorm) max(totalPhotonsNorm)]);

    hold off;
    box on
    set(gca, 'YLim', [oiWavelengthAxis(1) oiWavelengthAxis(end)], ...
        'XLim', [oiTimeAxis(1) oiTimeAxis(end)], 'FontSize', 14);
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('wavelength (nm)', 'FontSize', 14, 'FontWeight', 'bold');
    title('optical image photon rate (image center)', 'FontSize', 14);
    hC = colorbar('Location', 'NorthOutside', 'parent', timeSeriesTab);
    hC.FontSize =  14;
    hC.Label.String = 'photons/sec';
    axis 'xy'
    colormap(gray(1024));

    %% Plot the eye movement sequence (different colors for different OIs)
    subplot('Position', [0.25 0.07 0.22 0.89]);
    hold on;
    eyeMovementRange = [-100 100];

    plot(absorptionsTimeAxis, ...
        eyeMovementSequence(:, 1) * theConeMosaic.pigment.width * 1e6, ...
        '.', 'MarkerSize', 15, 'Color', 'r');
    hold on;
    plot(absorptionsTimeAxis, ...
        eyeMovementSequence(:, 2) * theConeMosaic.pigment.height * 1e6, ...
        '.', 'MarkerSize', 15, 'Color', 'b');
    plotStimulusTimes(eyeMovementRange);

    box on
    set(gca, 'YLim', [eyeMovementRange(1) eyeMovementRange(end)], ...
        'XLim', [oiTimeAxis(1) oiTimeAxis(end)], 'FontSize', 14);
    legend({'eye position (X)', 'eye position (Y)'});
    ylabel('X, Y eye position (microns)', 'FontSize', 14, ...
        'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('X, Y eye movements', 'FontSize', 14);

    %% Plot the LMS isomerizations
    if (theConeMosaic.rows ==1) && (theConeMosaic.cols == 3)
       referenceConeRows = [1 1 1]; referenceConeCols = [1 2 3];
    else
       % Find the (row, col) coords of the center-most L, M and S-cone
       for k = 1:3
            coneIndices = find(theConeMosaic.pattern == k + 1);
            [~, idx] = ...
                min(sum((theConeMosaic.coneLocs(coneIndices, :)) .^ 2, 2));
            [referenceConeRows(k), referenceConeCols(k)] = ...
                ind2sub(size(theConeMosaic.pattern), coneIndices(idx));
       end
    end
    subplot('Position', [0.50 0.07 0.22 0.89]);
    isomerizationRange = [min(isomerizationRateSequence(:)), ...
        1.05 * max(isomerizationRateSequence(:))];
    hold on
    coneColors = [1 0 0; 0 1 0; 0 0 1];
    for k = 1:3
        plot(absorptionsTimeAxis, squeeze(isomerizationRateSequence(1, ...
            referenceConeRows(k), referenceConeCols(k), :)), '.', ...
            'Color', squeeze(coneColors(k, :)), 'MarkerSize', 15, ...
            'LineWidth', 1.5);
    end
    plotStimulusTimes(isomerizationRange);

    hold off;
    set(gca, 'YLim', isomerizationRange, ...
        'XLim', [oiTimeAxis(1) oiTimeAxis(end)], 'FontSize', 14);
    ylabel('isomerization rate (R*/cone/sec)', 'FontSize', 14, ...
        'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('L, M, S-cone isomerization rates', 'FontSize', 14);

    %% Plot the photocurrents
    subplot('Position', [0.75 0.07 0.22 0.89]);
    photoCurrentRange = [min(photoCurrentSequence(:)), ...
        max(photoCurrentSequence(:)) + 2];
    hold on;
    for k = 1:3
        plot(responseTimeAxis, squeeze(photoCurrentSequence(1, ...
            referenceConeRows(k), referenceConeCols(k), :)), 'k.', ...
            'Color', squeeze(coneColors(k, :)), 'MarkerSize', 15, ...
            'LineWidth', 1.5);
    end
    plotStimulusTimes(photoCurrentRange);

    hold off;
    set(gca, 'XLim', [oiTimeAxis(1) oiTimeAxis(end)], ...
        'YLim', photoCurrentRange, 'FontSize', 14);
    ylabel('photocurrent (pA)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('time (seconds)', 'FontSize', 14, 'FontWeight', 'bold');
    title('@osLinear response', 'FontSize', 14);

    % Switch to eye movements tab
    set(tabGroup, 'SelectedTab', eyeMovementsTab);

    drawnow

    function plotStimulusTimes(signalRange)
    % Plot lines demarkating each OI time duration
    %
    % Syntax:
    %   plotStimulusTimes(signalRange)
    %
    % Description:
    %    A twice-embedded function inside plotEverything, designed to plot
    %    the lines demarcating each OI time duration.
    %
    % Inputs:
    %    signalRange -
        for oiIdx = 1:theOIsequence.length
            plot(oiTimeAxis(oiIdx)*[1 1], signalRange, 'k-');
        end
        % Plot the origin in magenta
        plot([0 0], signalRange, '-', 'Color', [0.7 0.1 0.3], 'LineWidth', 2);
    end

end

function condData = makeConditionSet(conditionSet)
% Helper function to create a condition set structure based on case set.
%
% Syntax:
%   condData = makeConditionSet(conditionSet)
%
% Description:
%    An embedded helper function in t_osTimeStep designed to create a
%    condition set structure out of a specified switch case set.
%
% Inputs:
%    conditionSet - Numeric. An integer for the switch case to determine
%                   which of the possible case structure sets to return.
%                   Condition sets explore varying effects. Options are:
%           1: Examine effects of varying integration time.
%           2: Examine effects of varying os.timeStep
%           3: Examine & contrast magnitudes of photon noise vs OS noise
%           4: Custom condition (customize!)
%
% Outputs:
%    condData     - Struct. A structure containing the conditions specified
%                   by the set case. Structure parameters include:
%       Steady: mosaicSize, meanLuminance, modulation, modulationRegion,
%               stimulusSamplingInterval, osTimeStep, integrationTime,
%               photonNoise, osNoise
%       Varied: (by case below)
%           1 - integrationTime
%           2 - osTimeStep
%           3 - photonNoise & osNoise
%           4 - stimulusSamplingInterval & integrationTime
%
% Optional key/value pairs:
%    None.
%

    % scene mean luminance
    meanLuminance = 1500;

    % Assemble conditions to run.
    condData = {};

    switch conditionSet
        case 1  % Effects of varying the integrationTime
            %% Steady parameters:
            % mosaicSize: NaN. Indicates a single L, M, and S cone.
            % meanLuminance: 1500. The scene mean luminance.
            % modulation: 0.5. A 50%  modulation against the background.
            % modulationRegion: 'Center'. Modulate the central image.
            %    Options are either Full or Center.
            % stimulusSamplingInterval: .02. A 50Hz stimulus refresh (i.e.
            %    20msec per optical image.)
            % osTimeStep: 0.1 millisecond. Time step per optical image.
            % integrationTime: NaN. Not yet set (see below).
            % photonNoise: 'Random'. Select photon noise from the
            %    following: {'random', 'frozen', 'none'}
            % osNoise: 'None'. Select OS noise from the following:
            %    {'random', 'frozen', 'none'}
            c0 = struct(...
                'mosaicSize', nan, ...
                'meanLuminance', meanLuminance, ...
                'modulation', 0.5, ...
                'modulationRegion', 'CENTER', ...
                'stimulusSamplingInterval', 1 / 50, ...
                'osTimeStep', 0.1 / 1000, ...
                'integrationTime', nan, ...
                'photonNoise', 'random', ...
                'osNoise', 'none');

            %% Varied params
            % 100ms is longer than stimulus sampling interval, therefore
            % < 1 eye movement / oi
            c0.integrationTime = 100 / 1000;
            condData{numel(condData) + 1} = c0;

            % 10 ms is smaller than the stimulus sampling interval,
            % therefore 2 eye movements / oi
            c0.integrationTime = 10 / 1000;
            condData{numel(condData) + 1} = c0;

            % 1 ms is also smaller, therefore 20 eye movements / oi
            c0.integrationTime = 1 / 1000;
            condData{numel(condData) + 1} = c0;

        case 2  % Examine effects of varying the os.timeStep
            %% Steady params
            % mosaicSize: NaN. Indicates a single L, M, and S cone.
            % meanLuminance: 1500. The scene mean luminance.
            % modulation: 0.5. A 50%  modulation against the background.
            % modulationRegion: 'Center'. Modulate the central image.
            %    Options are either Full or Center.
            % stimulusSamplingInterval: .1. A 100Hz stimulus refresh.
            % osTimeStep: NaN. Not yet set (see below).
            % integrationTime: .02. Integration time of 20 milliseconds.
            % photonNoise: 'Random'. Select photon noise from the
            %    following: {'random', 'frozen', 'none'}
            % osNoise: 'None'. Select OS noise from the following:
            %    {'random', 'frozen', 'none'}
            c0 = struct(...
                'mosaicSize', nan, ...
                'meanLuminance', meanLuminance, ...
                'modulation', 0.5, ...
                'modulationRegion', 'CENTER', ...
                'stimulusSamplingInterval', 1 / 10, ...
                'integrationTime', 20 / 1000, ...
                'osTimeStep', nan, ...
                'photonNoise', 'random', ...
                'osNoise', 'none');

            %% Varied params
            % 1 milliseconds.
            c0.osTimeStep = 1 / 1000;
            condData{numel(condData) + 1} = c0;

            % 0.5 milliseconds.
            c0.osTimeStep = 0.5 / 1000;
            condData{numel(condData) + 1} = c0;

            % 0.1 milliseconds.
            c0.osTimeStep = 0.1 / 1000;
            condData{numel(condData) + 1} = c0;

        case 3  % Examine & contrast magnitudes of photon noise vs OS noise

            %% Steady params
            % mosaicSize: NaN. Indicates a single L, M, and S cone.
            % meanLuminance: 1500. The scene mean luminance.
            % modulation: 0.5. A 50%  modulation against the background.
            % modulationRegion: 'Center'. Modulate the central image.
            %    Options are either Full or Center.
            % stimulusSamplingInterval: .04. A 25Hz stimulus refresh.
            % osTimeStep: 0.1 millisecond. Time step per optical image.
            % integrationTime: 5 / 1000. Integration time of 5 milliseconds
            % photonNoise: NaN. Not yet set (see below).
            % osNoise: NaN. Not yet set (see below).
            c0 = struct(...
                'mosaicSize', nan, ...
                'meanLuminance', meanLuminance, ...
                'modulation', 0.5, ...
                'modulationRegion', 'CENTER', ...
                'stimulusSamplingInterval', 1 / 10, ...
                'integrationTime', 20 / 1000, ...
                'osTimeStep', 0.1 / 1000, ...
                'photonNoise', nan, ...
                'osNoise', nan);

            %% Varied params
            % No noise (Options are: 'random', 'frozen', 'none')
            c0.photonNoise = 'none';
            c0.osNoise = 'none';
            condData{numel(condData) + 1} = c0;

            % Photon noise (Options are: 'random', 'frozen', 'none')
            c0.photonNoise = 'random';
            c0.osNoise = 'none';
            condData{numel(condData) + 1} = c0;

            % OS noise (Options for either are: 'random', 'frozen', 'none')
            c0.photonNoise = 'none';
            c0.osNoise = 'random';
            condData{numel(condData) + 1} = c0;

            % Photon & OS noise (Options are: 'random', 'frozen', 'none')
            c0.photonNoise = 'random';
            c0.osNoise = 'random';
            condData{numel(condData) + 1} = c0;

        case 4  % Custom condition - customize it!
            % Customizations
            stimRefreshRateInHz = 66;
            integrationTime = 50 / 1000;

            %% Steady parameters
            % mosaicSize: 0.5. A Field of view of 0.5 x 0.5 degrees.
            % meanLuminance: 1500. The scene mean luminance.
            % modulation: 0.5. A 50%  modulation against the background.
            % modulationRegion: 'Center'. Modulate the central image.
            %    Options are either Full or Center.
            % stimulusSamplingInterval: 1 / 66. A 66Hz stimulus refresh.
            % osTimeStep: 0.1 millisecond. Time step per optical image.
            % integrationTime: .05. Integration time of 50 milliseconds
            % photonNoise: 'Random'. Select photon noise from the
            %    following: {'random', 'frozen', 'none'}
            % osNoise: 'None'. Select OS noise from the following:
            %    {'random', 'frozen', 'none'}
            c0 = struct(...
                'mosaicSize', 0.5, ...
                'meanLuminance', meanLuminance, ...
                'modulation', 0.5, ...
                'modulationRegion', 'CENTER', ...
                'stimulusSamplingInterval', 1 / stimRefreshRateInHz, ...
                'integrationTime', integrationTime, ...
                'osTimeStep', 0.1/1000, ...
                'photonNoise', 'none', ...
                'osNoise', 'none');

           condData{numel(condData) + 1} = c0;
    end
end
