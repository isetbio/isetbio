function t_conesMosaicAbsorptionOISequenceEcc
% Demonstrate applying ecc-based efficiency correction for dynamic stimulus
%
% Syntax:
%   t_conesMosaicAbsorptionOISequenceEcc
%
% Description:
%    Demonstrates the effect of applying eccentricity-based cone efficiency
%    correction for a dynamic (ramped on/off) stimulus.
%
%    This demo includes the default fixational eye movements which may
%    contain microsaccades.
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
%    06/16/18  NPC  ISETBIO Team, 2018
%    06/16/18  npc  Wrote it.
%    08/30/18  jnm  Formatting

%% Initialize
ieInit;

%% Set an optical image sequence
% The sequence corresponds to a Gabor stimulus whose contrast is ramped up
% and down using a Gaussian modulation function.
stimFOV = 0.8;
stimDurationSeconds = 0.5;
sigmaSeconds = 0.15;
[theOIsequence, theScene, thePeakOI] = ...
    generateTheOISequence(stimFOV, stimDurationSeconds, sigmaSeconds);
peakPhotons = [];

%% Location of mosaic to be used
mosaicDir = ...
    strrep(isetRootPath, 'isettools', 'tutorials/t_cones/resources');
mosaicFileName = fullfile(mosaicDir, 'mosaic.mat');

%% Build a cone mosaic or load existing one
loadPreviouslySavedMosaic = true;
if (~loadPreviouslySavedMosaic)
    mosaicFOV = 0.6;
    cMosaic  = coneMosaicHex(9, ...
            'name', 'test', ...
            'fovDegs', mosaicFOV, ...
            'eccBasedConeDensity', true, ...
            'sConeMinDistanceFactor', 2.0, ...
            'sConeFreeRadiusMicrons', 0.15 * 300, ...
            'spatialDensity', [0 6/10 3/10 1/10], ...
            'latticeAdjustmentPositionalToleranceF', 0.5, ...
            'latticeAdjustmentDelaunayToleranceF', 0.05, ...
            'maxGridAdjustmentIterations', 200, ...
            'marginF', []);
    save(mosaicFileName, 'cMosaic', '-v7.3');
else
    load(mosaicFileName, 'cMosaic');
    mosaicFOV = max(cMosaic.fov);
end

% No photon noise because we want to see the mean isomerizations
cMosaic.noiseFlag = 'none';

% 1 msec integration time
cMosaic.integrationTime = 2 / 1000;

% No correction for ecc-dependent change in cone efficiency
cMosaic.eccBasedConeQuantalEfficiency = false;

%% Compute the number of eye movements
% This calculation is based on this mosaic's integration time and the
% current oiSequence.
eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(...
    cMosaic.integrationTime);

% Generate eye eyemovements for 1 trial
nTrials = 1;
emPaths = cMosaic.emGenSequence(eyeMovementsNum, 'nTrials', nTrials);

% Save the eye movement paths so we can apply it to second mosaic
theEmPositions = emPaths;

% Make an identical copy of first mosaic
cMosaic2 = cMosaic.copy();

% but with ecc-based cone efficiency corrections
cMosaic2.eccBasedConeQuantalEfficiency = true;

%% Compute absorptions for both mosaics
absorptions1 = cMosaic.computeForOISequence(theOIsequence, ...
        'currentFlag', false, 'emPaths', theEmPositions);

absorptions2 = cMosaic2.computeForOISequence(theOIsequence, ...
        'currentFlag', false, 'emPaths', theEmPositions);

% Extract response time axis
timeAxisMsec = 1000 * cMosaic.timeAxis;

% Extract the optical image time axis
oiTimeAxisMsec = 1000 * theOIsequence.timeAxis();

% Extract absorptions range across both mosaics
signalRange = [...
    min([min(absorptions1(:)) min(absorptions2(:))])...
    max([max(absorptions1(:)) max(absorptions2(:))])];

%% Initialize figure and video
generateVideo = true;
[hFig, subplotPosVectors, writerObj] = initFigure(generateVideo);

% Set the visualizedFOV
visualizedFOV = mosaicFOV;

%% Plot the two mosaics
plotMosaic(cMosaic, subplotPosVectors, [1 1], visualizedFOV, ...
    'mosaic with constant cone efficiency');
plotMosaic(cMosaic2, subplotPosVectors, [2 1], visualizedFOV, ...
    'mosaic with ecc-based cone efficiency');

%% Generate LUTs to represent mosaic activation
activationLUT = brewermap(1024, '*greys');
diffActivationLUT = brewermap(1024, '*RdBu');

%% Plot each time bin of the two mosaic activations (and their difference)
for visualizedTrialIndex = 1:size(theEmPositions, 1)
    for timeBin = 1:size(theEmPositions, 2)
        %% Extract the activation patterns for current trial and time bin
        visActivationPattern1 = ...
            squeeze(absorptions1(visualizedTrialIndex, :, timeBin));
        visActivationPattern2 = ...
            squeeze(absorptions2(visualizedTrialIndex, :, timeBin));
        %% Reshape into full 2D pattern
        visActivationPattern1 = ...
            cMosaic.reshapeHex1DmapToHex2Dmap(visActivationPattern1);
        visActivationPattern2 = ...
            cMosaic.reshapeHex1DmapToHex2Dmap(visActivationPattern2);
        differentialActivationPattern = ...
            visActivationPattern2 - visActivationPattern1;

        %% Plot activation of first mosaic
        plotMosaicActivation(cMosaic, visActivationPattern1, ...
            signalRange, theEmPositions, visualizedTrialIndex, timeBin, ...
            activationLUT, [0 0 0], visualizedFOV, subplotPosVectors, ...
            [1 2], 'absorptions in constant cone efficiency mosaic', ...
            'displaylabelX', false, 'displaylabelY', false, ...
            'displayTicksX', false, 'displayTicksY', false);

        %% Plot activation of second mosaic
        plotMosaicActivation(cMosaic2, visActivationPattern2, ...
            signalRange, theEmPositions, visualizedTrialIndex, timeBin, ...
            activationLUT, [0 0 0], visualizedFOV, subplotPosVectors, ...
            [2 2], 'absorptions in ecc-based cone efficiency mosaic', ...
            'displaylabelY', false, 'displayTicksY', false);

        %% Plot the difference in absorptions b/n the 2 mosaics
        plotMosaicActivation(cMosaic2, differentialActivationPattern, ...
            round(max(signalRange) / 1.8) * [-1 1], theEmPositions, ...
            visualizedTrialIndex, timeBin, diffActivationLUT, [1 1 1], ...
            visualizedFOV, subplotPosVectors, [1 3], ...
            sprintf('absorptions difference (trial:%d, t:%2.1f ms)', ...
            visualizedTrialIndex, timeAxisMsec(timeBin)), ...
            'displaylabelX', false, 'displayTicksX', false, ...
            'displaylabelY', false, 'displayTicksY', false);

        %% Plot the current optical image
        [~, idx] = min(abs(oiTimeAxisMsec - timeAxisMsec(timeBin)));
        theCurrentOI = theOIsequence.frameAtIndex(idx);
        peakPhotons = plotOpticalImage(theCurrentOI, thePeakOI, ...
            subplotPosVectors, [2 3], visualizedFOV, ...
            timeAxisMsec(timeBin), peakPhotons);

        %% Render frame
        drawnow;
        if (generateVideo)
            % Add frame to video stream
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end

    end % timeBin
end % visualizedTrialIndex

% Close video stream
if (generateVideo), close(writerObj); end

end

% ---------------- SUPPORT FUNCTIONS --------------------
function [hFig, subplotPosVectors, writerObj] = initFigure(generateVideo)
% initialize figure
%
% Syntax:
%   [hFig, subplotPosVectors, writerObj] = initFigure(generateVideo)
%
% Description:
%    A support function for t_conesMosaicAbsorptionOISequenceEcc that
%    initializes the figure.
%
% Inputs:
%    generateVideo     - Boolean. A boolean indicating whether or not to
%                        generate the video.
%
% Outputs:
%    hFig              - Handle. The handle for the generated figure.
%    subplotPosVectors - Struct. A structure containing the subplot
%                        position vectors.
%    writerObj         - Object. A video writer object.
%
% Optional key/value pairs:
%    None.
%
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 3, ...
       'heightMargin', 0.03, ...
       'widthMargin', 0.03, ...
       'leftMargin', 0.01, ...
       'rightMargin', 0.03, ...
       'bottomMargin', 0.04, ...
       'topMargin', 0.01);

    % Generate figure
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1295 792], 'Color', [1 1 1]);

    if (generateVideo)
        videoFile = fullfile(tempdir, 'video.m4v');
        writerObj = VideoWriter(videoFile, 'MPEG-4'); % H264 format
        writerObj.FrameRate = 30;
        writerObj.Quality = 100;
        % Open video stream
        open(writerObj);
    else
        writerObj = [];
    end
end

function plotScene(scene, subplotPosVectors, subplotPos, visualizedFOV)
% plot the scene
%
% Syntax:
%   plotScene(scene, subplotPosVectors, subplotPos, visualizedFOV)
%
% Description:
%    A support function for t_conesMosaicAbsorptionOISequenceEcc that
%    allows you to plot the scene.
%
% Inputs:
%    scene             - Struct. A scene structure.
%    subplotPosVectors - Struct. A structure containing the subplot
%                        position vectors.
%    subplotPos        - Vector. The subplot position vector. This is the
%                        (x,y) position of the subplot within the figure.
%    visualizedFOV     - Numeric. The visualized field of view.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    stimFOV =  sceneGet(scene, 'fov');
    % Plot the scene
    ax = subplot('Position', ...
        subplotPosVectors(subplotPos(1), subplotPos(2)).v);
    xspace = (1:sceneGet(scene, 'cols'));
    xspace = xspace - mean(xspace);
    xspace = xspace / max(xspace) * stimFOV / 2;
    yspace = (1:sceneGet(scene, 'rows'));
    yspace = yspace - mean(yspace);
    yspace = yspace / max(yspace) * stimFOV / 2;
    luminanceMap = sceneGet(scene, 'luminance');

    maxLuminance = ceil(max(luminanceMap(:)));
    minLuminance = 0;

    imagesc(ax, xspace, yspace, ...
        (luminanceMap - minLuminance) / (maxLuminance-minLuminance));
    set(gca, 'CLim', [0 1]);
    axis('image');
    hold on;
    plot([xspace(1) xspace(end)], [0 0], 'g-', 'LineWidth', 1.5);
    plot([0 0], [yspace(1) yspace(end)], 'g-', 'LineWidth', 1.5);
    hold off;

    xlabel(ax, 'space (degs)');
    lut = gray(1024);
    colormap(ax, lut);
    set(ax, 'XTick', visualizedFOV / 2 * [-1 0 1], ...
        'YTick', visualizedFOV / 2 * [-1 0 1], ...
        'XLim', visualizedFOV / 2 * [-1 1], ...
        'YLim', visualizedFOV / 2 * [-1 1], 'YTickLabels', {}, ...
        'FontSize', 14, 'Color', [0 0 0]);
    xtickformat(ax, '%0.1f');

    lumTicks = 0:0.1:1.0;
    hC = colorbar();
    hC.Ticks = lumTicks;
    hC.TickLabels = sprintf('%2.1f\n', ...
        (lumTicks * (maxLuminance - minLuminance) + minLuminance));
    hC.Label.String = 'luminance (cd/m2)';
    title(ax, 'scene luminance');
end

function peakPhotons = plotOpticalImage(oi, thePeakOI, ...
    subplotPosVectors, subplotPos, visualizedFOV, currentTime, peakPhotons)
% plot the optical image
%
% Syntax:
%   peakPhotons = plotOpticalImage(oi, thePeakOI, subplotPosVectors, ...
%       subplotPos, visualizedFOV, currentTime, peakPhotons)
%
% Description:
%    A support function for t_conesMosaicAbsorptionOISequenceEcc that
%    allows you to plot the optical image.
%
% Inputs:
%    oi                - Struct. An optical image structure.
%    thePeakOI         - Struct. An optical image structure.
%    subplotPosVectors - Struct. The structure containing the subplot
%                        position vectors.
%    subplotPos        - Vector. The subplot position vector. This is the
%                        (x,y) position of the subplot within the figure.
%    visualizedFOV     - Numeric. The visualized field of view.
%    currentTime       - Numeric. The current time bin.
%    peakPhotons       - Numeric. The number of peak photons.
%
% Outputs:
%    peakPhotons      - Numeric. The number of peak photons.
%
% Optional key/value pairs:
%    None.
%
    stimFOV =  oiGet(oi, 'fov');
    % Plot the scene
    ax = subplot('Position', ...
        subplotPosVectors(subplotPos(1), subplotPos(2)).v);
    xspace = (1:oiGet(oi, 'cols'));
    xspace = xspace - mean(xspace);
    xspace = xspace / max(xspace) * stimFOV / 2;
    yspace = (1:oiGet(oi, 'rows'));
    yspace = yspace - mean(yspace);
    yspace = yspace / max(yspace) * stimFOV / 2;

    % Get the retinal photons map
    retinalPhotonsMap = oiGet(oi, 'photons');

    % Get the retinal photons map at 550 nm
    wave = oiGet(oi, 'wave');
    targetWave = 550;
    [~, idx] = min(abs(wave - targetWave));
    targetWave = wave(idx);
    retinalPhotonsMap550nm = squeeze(retinalPhotonsMap(:, :, idx));

    if (isempty(peakPhotons))
        retinalPeakPhotonsMap = oiGet(thePeakOI, 'photons');
        retinalPeakPhotonsMap = squeeze(retinalPeakPhotonsMap(:, :, idx));
        peakPhotons = ceil(max(retinalPeakPhotonsMap(:)));
    end

    maxPhotons = peakPhotons;
    minPhotons = 0;

    imagesc(ax, xspace, yspace, ...
        (retinalPhotonsMap550nm - minPhotons) / (maxPhotons - minPhotons));
    set(ax, 'CLim', [0 1]);
    axis('image');
    hold on;
    plot([xspace(1) xspace(end)], [0 0], 'g-', 'LineWidth', 1.5);
    plot([0 0], [yspace(1) yspace(end)], 'g-', 'LineWidth', 1.5);
    hold off;

    xlabel(ax, 'space (degs)');

    lut = gray(1024);
    colormap(gca, lut);
    set(ax, 'XTick', visualizedFOV / 2 * [-1 0 1], ...
        'YTick', visualizedFOV / 2 * [-1 0 1], ...
        'XLim', visualizedFOV / 2 * [-1 1], ...
        'YLim', visualizedFOV / 2 * [-1 1], 'YTickLabels', {}, ...
        'FontSize', 14, 'Color', [0 0 0]);
    xtickformat(ax, '%0.1f');

    photonsTicks = 0:0.1:1.0;
    scaleFactor = 1e9;
    scaleFactorSting = '1E9';
    hC = colorbar();
    hC.Ticks = photonsTicks;
    hC.TickLabels = sprintf('%2.1f\n', (photonsTicks / 1e3 * ...
        (maxPhotons - minPhotons) + minPhotons) / scaleFactor);
    hC.Label.String = ...
        sprintf('photon rate (photons / msec) x %s', scaleFactorSting);

    title(ax, sprintf('retinal photon rate @%2.0f nm (%2.0f msec)', ...
        targetWave, currentTime));
end

function plotMosaic(cMosaic, subplotPosVectors, subplotPos, ...
    visualizedFOV, titleText)
% Plot the mosaic
%
% Syntax:
%   plotMosaic(cMosaic, subplotPosVectors, subplotPos, ...
%       visualizedFOV, titleText)
%
% Description:
%    A support function for t_conesMosaicAbsorptionOISequenceEcc to plot
%    the cone mosaic.
%
% Inputs:
%    cMosaic           - Object. A cone mosaic object.
%    subplotPosVectors - Struct. A structure containing the subplot
%                        position vectors.
%    subplotPos        - Vector. The subplot position vector. This is the
%                        (x,y) position of the subplot within the figure.
%    visualizedFOV     - Numeric. The visualized field of view.
%    titleText         - String. The plot title string.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
    micronPerDeg = cMosaic.micronsPerDegree;
    ax = subplot('Position', ...
        subplotPosVectors(subplotPos(1), subplotPos(2)).v);
    cMosaic.visualizeGrid(...
        'axesHandle', ax, ...
        'visualizedConeAperture', 'lightCollectingArea', ...
        'apertureShape', 'disks', ...
        'overlayHexMesh', true);
    set(ax, 'XLim', visualizedFOV / 2 * [-1 1] * micronPerDeg * 1e-6, ...
         'YLim', visualizedFOV / 2 * [-1 1] * micronPerDeg * 1e-6, ...
         'XTick', [], 'YTick', [], 'FontSize', 14);
    title(ax, titleText);
    colormap(ax, gray(1024));
    hC = colorbar(ax);
    hC.Ticks = [];
end

function plotMosaicActivation(cMosaic, visActivationPattern, ...
    signalRange, theEmPositions, visualizedTrialIndex, timeBin, ...
    activationLUT, backgroundColor, visualizedFOV, subplotPosVectors, ...
    subplotPos, titleText, varargin)
% plot the mosaic activation
%
% Syntax:
%   plotMosaicActivation(cMosaic, visActivationPattern, ...
%       signalRange, theEmPositions, visualizedTrialIndex, timeBin, ...
%       activationLUT, backgroundColor, visualizedFOV, ...
%       subplotPosVectors, subplotPos, titleText)
%
% Description:
%    A support function for t_conesMosaicAbsorptionOISequenceEcc that
%    allows you to plot the mosaic activation.
%
% Inputs:
%    cMosaic              - Object. A cone mosaic object.
%    visActivationPattern - Matrix. The visualized activation pattern.
%    signalRange          - Vector. A 1x2 vector containing signal range.
%    theEmPositions       - Matrix. A matrix containing the eye movement
%                           positions for the R & L eyes.
%    visualizedTrialIndex - Numeric. The trial index.
%    timeBin              - Numeric. Which time bin you are in.
%    activationLUT        - Matrix. The activation LUT.
%    backgroundColor      - Vector. The [R, G, B] numeric vector
%                           representing the subplot background color.
%    visualizedFOV        - Numeric. The visualized field of view.
%    subplotPosVectors    - Struct. The subplot position vectors.
%    subplotPos           - Vector. A 1x2 vector containing the
%                           subplot positions.
%    titleText            - String. The Figure title as a string.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    displayLabelX        - Boolean. Whether or not to display the X axis
%                           label. Default true.
%    displayLabelY        - Boolean. Whether or not to display the Y axis
%                           label. Default true.
%    displayTicksX        - Boolean. Whether or not to display the X axis
%                           tick marks. Default true.
%    displayTicksY        - Boolean. Whether or not to display the Y axis
%                           tick marks. Default true.
%
    p = inputParser;
    p.addParameter('displaylabelX', true, @islogical);
    p.addParameter('displaylabelY', true, @islogical);
    p.addParameter('displayTicksX', true, @islogical);
    p.addParameter('displayTicksY', true, @islogical);
    p.parse(varargin{:});

    ax = subplot('Position', ...
        subplotPosVectors(subplotPos(1), subplotPos(2)).v);
    cMosaic.renderActivationMap(ax, visActivationPattern, ...
             'signalRange', signalRange, ...
             'visualizedConeAperture', 'geometricArea', ...
             'mapType', 'modulated disks', ...
             'showColorBar', true, ...
             'labelColorBarTicks', true, ...
             'titleForColorBar', sprintf('R* / cone /%2.0fmsec', ...
             1000 * cMosaic.integrationTime), ...
             'colorMap', activationLUT, ...
             'crossHairPosition', squeeze(theEmPositions(...
             visualizedTrialIndex, timeBin, :)) * ...
             cMosaic.patternSampleSize(1), ...
             'visualizedFOV', visualizedFOV, ...
             'backgroundColor', backgroundColor);
    set(ax, 'FontSize', 14);
    if (~p.Results.displaylabelX), xlabel(ax, ''); end
    if (~p.Results.displaylabelY), ylabel(ax, ''); end
    if (~p.Results.displayTicksX), set(ax, 'XTickLabels', {}); end
    if (~p.Results.displayTicksY), set(ax, 'YTickLabels', {}); end
    title(ax, titleText);
end

function [theOIsequence, theScene, thePeakOI] = ...
    generateTheOISequence(stimFOV, stimDurationSeconds, sigmaSeconds)
% Generate the OI sequence
%
% Syntax:
%   [theOIsequence, theScene, thePeakOI] = ...
%       generateTheOISequence(stimFOV, stimDurationSeconds, sigmaSeconds)
%
% Description:
%    A support function for t_conesMosaicAbsorptionOISequenceEcc that
%    allows you to generate the optical image sequence.
%
% Inputs:
%    stimFOV             - Numeric. The stimulus field of view.
%    stimDurationSeconds - Numeric. The stimulus duration in seconds.
%    sigmaSeconds        - Numeric. The standard deviation of seconds.
%
% Outputs:
%    theOISequence       - Object. The optical image sequence object.
%    theScene            - Struct. The associated scene structure.
%    thePeakOI           - Struct. The peak optical image structure.
%
% Optional key/value pairs:
%    None.
%
    sparams.fov = stimFOV;    % scene field of view in degrees
    sparams.distance = 0.57;  % viewing distance in meters

    % Gabor parameters
    P(2) = harmonicP;
    P(2).ph = 0;
    P(2).freq = 8;
    P(2).row = 128;
    P(2).col = 128;
    P(2).GaborFlag = .12;
    P(1) = P(2);
    P(1).contrast = 0;

    % Stimulus temporal properties
    stimRefreshRate = 10 / 1000;  % 100 Hz display
    stimSampleTimes = (0:ceil(stimDurationSeconds / stimRefreshRate)) ...
        * stimRefreshRate;
    sigma = sigmaSeconds / stimRefreshRate;
    stimModulationEnvelope = ieScale(fspecial('gaussian', ...
        [1, numel(stimSampleTimes)], sigma), 0, 1);

    % Make an optical image sequence of the stimulus
    [theOIsequence, theScene] = oisCreate('harmonic', ...
        'blend', stimModulationEnvelope, ...
        'sampleTimes', stimSampleTimes, ...
        'testParameters', P, 'sceneParameters', sparams);

    % Get the oi at the peak of the modulation
    [~, idx] = max(stimModulationEnvelope);
    thePeakOI = theOIsequence.frameAtIndex(idx);
end

