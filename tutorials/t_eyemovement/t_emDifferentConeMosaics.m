function t_emDifferentConeMosaics
%%t_emDifferentConeMosaics  Tests whether the eye movement dynamics are independent of the spatial characteristics of the cone mosaic.  
%
% Description:
%    Demonstrates eye movements dynamics for three types of cone mosaics: 
%     (a) a rectangular cone mosaic,
%     (b) a hexagonal cone mosaic with regular cone spacing and larger (than the default values in isetbio) aperture & spacing 
%     (c) a hexagonal cone mosaic with eccentricity-based cone spacing and smaller than the default value cone aperture
%
%    Each of the three rows in the generated figure corresponds to one of the three cone mosaics studied.
%    The first column  plots the cone mosaic with the eye movement from the first trial superimposed in red.
%    The second column plots the x- and y-components of the eye movement from the first trial.
%    The third column plots the mean spectra of the x- and y-components across all computed trials (256 here)
%
%    Notes: The eye movements appear to have similar dynamics across all mosaics. 
%      They have similar x- and y-ranges, and their temporal spectra have
%      similar shape except for a vertical shift: the power for the rect-mosaic is lower. 
%      This is because the em positions are quantized to the mosaic pattern sample size and the sample size in
%      the hex mosaic is typically set to 4-9 times smaller than in the rect mosaic (resamplingFactor).
%      Indeed, if we set the resamplingFactor to 1, the spectra of the rect and hex mosaics agree pretty well.
%
%    Overall, this demonstrates that the eye movement dynamics are independent of the employed mosaic's 
%    spatial characteristics (pattern size, the cone spacing, and cone aperture).
%
%
% NPC, ISETBIO Team, 2017
%
% 10/03/17  npc  Wrote it.
%
%% Initialize
ieInit;

%% Set the random seed
rng(1);

fovDegs = 0.15;
integrationTimeMillisecs = 1;
nTrials = 256;
eyeMovementsPerTrial = 10100;

%% The default rect mosaic
[cm, theEMpaths, emSpectrumXo, emSpectrumYo, tfAxis] = t_emRectMosaic(fovDegs, integrationTimeMillisecs, nTrials, eyeMovementsPerTrial);
% Plot the eye movements during the first trial
iTrial = 1;
plotEMs(1, cm, theEMpaths, iTrial, emSpectrumXo, emSpectrumYo, [], [], tfAxis);

%% A regular hex cone mosaic with large separation (customLambda: 7 microns) and large inner segment diameter (4 microns)
resamplingFactor = 6;
eccBasedConeDensity = false; customLamda = 7; customInnerSegmentDiameter = 4;
[cm, theEMpaths, emSpectrumX, emSpectrumY, tfAxis] = t_emHexMosaic(fovDegs, integrationTimeMillisecs, eccBasedConeDensity, customLamda, customInnerSegmentDiameter, resamplingFactor, nTrials, eyeMovementsPerTrial);
% Plot the eye movements during the first trial
iTrial = 1;
plotEMs(2, cm, theEMpaths, iTrial, emSpectrumX, emSpectrumY, emSpectrumXo, emSpectrumYo, tfAxis);

%% An ecc-based hex cone mosaic with small inner segment diameter (1 micron)
eccBasedConeDensity = true; customInnerSegmentDiameter = 1; customLamda = [];
[cm, theEMpaths, emSpectrumX, emSpectrumY, tfAxis] = t_emHexMosaic(fovDegs, integrationTimeMillisecs, eccBasedConeDensity, customLamda, customInnerSegmentDiameter, resamplingFactor, nTrials, eyeMovementsPerTrial);
% Plot the eye movements during the first trial
iTrial = 1;
plotEMs(3, cm, theEMpaths, iTrial, emSpectrumX, emSpectrumY, emSpectrumXo, emSpectrumYo, tfAxis);
end

function [cm, theEMpaths, emSpectrumX, emSpectrumY, tfAxis] = t_emRectMosaic(fovDegs, integrationTimeMillisecs, nTrials, eyeMovementsPerTrial)
% Generate the rect mosaic
cm = coneMosaic();
cm.setSizeToFOV(fovDegs);            

% Set the integration time
cm.integrationTime = integrationTimeMillisecs/1000;  

% Generate eye movement paths
theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
for iTrial= 1:nTrials
    theEMpaths(iTrial, :,:) = cm.emGenSequence(eyeMovementsPerTrial);
end

% Compute spectra
[emSpectrumX, emSpectrumY, tfAxis] = computeEMspectrum(cm.timeAxis*1000, theEMpaths);
end

function [cm, theEMpaths, emSpectrumX, emSpectrumY, tfAxis] = t_emHexMosaic(fovDegs, integrationTimeMillisecs, eccBasedConeDensity, customLamda, customInnerSegmentDiameter, resamplingFactor, nTrials, eyeMovementsPerTrial)
% Generate the hex mosaic
cm = coneMosaicHex(resamplingFactor, ...
        'fovDegs', fovDegs, ...
        'eccBasedConeDensity', eccBasedConeDensity, ... 
        'customLambda', customLamda, ...
        'customInnerSegmentDiameter', customInnerSegmentDiameter, ...
        'latticeAdjustmentPositionalToleranceF', 0.1, ...
        'latticeAdjustmentDelaunayToleranceF', 0.1 ...
    );
% Set the integration time
cm.integrationTime = integrationTimeMillisecs/1000; 

% Generate eye movement paths
theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
for iTrial= 1:nTrials
    theEMpaths(iTrial, :,:) = cm.emGenSequence(eyeMovementsPerTrial);
end

% Compute spectra
[emSpectrumX, emSpectrumY, tfAxis] = computeEMspectrum(cm.timeAxis*1000, theEMpaths);
end


% ----- Spectrum analysis routine ----
function [emSpectrumX, emSpectrumY, tfAxis] = computeEMspectrum(timeAxisMillisecs, theEMpaths)
iCounter = 0;
for emTrial = 1:size(theEMpaths,1)
    iCounter = iCounter + 1;
    emSpectrumX(iCounter,:) = abs(fft(squeeze(theEMpaths(emTrial,:,1))));
    emSpectrumY(iCounter,:) = abs(fft(squeeze(theEMpaths(emTrial,:,2))));
end
emSpectrumX = mean(emSpectrumX, 1);
%x_ftspectrum = x_ftspectrum / max(x_ftspectrum);
emSpectrumY = mean(emSpectrumY, 1);
%y_ftspectrum = y_ftspectrum / max(y_ftspectrum);
maxF = 1000/(2*(timeAxisMillisecs(2)-timeAxisMillisecs(1)));
deltaF = maxF / (0.5*numel(timeAxisMillisecs));
tfAxis = (1:size(theEMpaths,2))*deltaF;
end


%----- PLOTTING ROUTINE -------
function plotEMs(figRow, cm, theEMpaths, iTrial, emSpectrumX, emSpectrumY, emSpectrumXo, emSpectrumYo, tfAxis)

% convert EMPath from cone index to microns
patternSampleSizeMicrons = cm.patternSampleSize(1)*1e6;
theEMpathsMicrons = theEMpaths * patternSampleSizeMicrons;
xRange = max(max(max(theEMpathsMicrons(:,:,1)))) - min(min(min(theEMpathsMicrons(:,:,1))));
yRange = max(max(max(theEMpathsMicrons(:,:,2)))) - min(min(min(theEMpathsMicrons(:,:,2))));

% convert cone positions to microns
if (isa(cm, 'coneMosaicHex'))
    conePositionsMicrons = cm.coneLocs*1e6/cm.resamplingFactor;
else
    conePositionsMicrons = cm.coneLocs*1e6;
end

% convert aperture size to microns
apertureSizeMicrons = cm.pigment.width*1e6;
spatialExtent = 24*[-1 1];
xTicks = (spatialExtent(1):4:spatialExtent(2));
yTicks = (spatialExtent(1):4:spatialExtent(2));
xTickLabels = sprintf('%2.0f\n', xTicks);
yTickLabels = sprintf('%2.0f\n', yTicks);
    
xoutline = cos((0:30:360)/180*pi) * apertureSizeMicrons/2;
youtline = sin((0:30:360)/180*pi) * apertureSizeMicrons/2;

% convert time to msec
timeAxisMillisecs = cm.timeAxis*1000;

if (figRow == 1)
    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1500 1300], 'Color', [1 1 1]);
end

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', 3, ...
           'heightMargin',  0.08, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.04, ...
           'rightMargin',    0.0, ...
           'bottomMargin',   0.06, ...
           'topMargin',      0.03);

% Plot eye movements on top of cone mosaic
subplot('Position', subplotPosVectors(figRow,1).v);
    
if (isa(cm, 'coneMosaicHex'))
    % Use the visualizeGrid() method to show the mosaic and the eye movemtns
    cm.visualizeGrid(...
        'axesHandle', gca, ...
        'overlayEMpath', squeeze(theEMpaths(1,:,:)), ...
        'overlayNullSensors', true, ...
        'apertureShape', 'disks', ...
        'visualizedConeAperture', 'both', ...
        'labelConeTypes', false,...
        'generateNewFigure', false);
    set(gca, 'XTick', xTicks*1e-6, 'YTick', yTicks*1e-6, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
    set(gca, 'XLim', spatialExtent*1e-6, 'YLim', spatialExtent*1e-6);
    customLambda = cm.customLambda;
    if (isempty(customLambda))
        customLambda = nan;
    end
    title(sprintf('Hex mosaic,  pattern sample size: %2.1f um\ncone sep: %g um, aperture: %2.1f um', cm.patternSampleSize(1)*1e6, customLambda, cm.pigment.pdWidth*1e6));
else
    hold on
    for k = 1:size(conePositionsMicrons,1)
        if (cm.pattern(k) > 1)
            xx = conePositionsMicrons(k,1) + xoutline;
            yy = conePositionsMicrons(k,2) + youtline;
            plot(xx,yy, 'k-');
        end
    end
    plot(theEMpathsMicrons(iTrial,:,1), theEMpathsMicrons(iTrial,:,2), 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', [1 0.5 0.5]);
    set(gca, 'XTick', xTicks, 'YTick', yTicks, 'XTickLabel', xTickLabels, 'YTickLabel', yTickLabels);
    set(gca, 'XLim', spatialExtent, 'YLim', spatialExtent);
    axis 'square'
    title(sprintf('Rect mosaic, pattern sample size: %2.1f um\ncone sep: %g um, aperture: %2.1f um', cm.patternSampleSize(1)*1e6, cm.patternSampleSize(1)*1e6, cm.pigment.pdWidth*1e6));
end
ylabel('space (microns)',  'FontWeight', 'bold');
if (figRow == 3)
    xlabel('space (microns)',  'FontWeight', 'bold');
end
set(gca, 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
box(gca, 'on'); grid(gca, 'off');
    

% Plot the time course of the x and y eye movement
subplot('Position', subplotPosVectors(figRow,2).v);
plot(timeAxisMillisecs, theEMpathsMicrons(iTrial,:,1), 'k-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerEdgeColor', 'm', 'Color', 'm');
hold on;
plot(timeAxisMillisecs, theEMpathsMicrons(iTrial,:,2), 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'Color', 'b');
set(gca, 'XLim', [timeAxisMillisecs(1) timeAxisMillisecs(end)], 'YLim', spatialExtent, 'YTick', yTicks);
grid on; box off;
set(gca, 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
hL = legend({'x-pos', 'y-pos'}, 'Location', 'SouthWest');
set(hL, 'FontSize', 14);
ylabel('space (microns)',  'FontWeight', 'bold'); 
if (figRow == 3)
    xlabel('time (msec)', 'FontWeight', 'bold'); 
end
title(sprintf('X/Y eye position sequence (trial #1)\nXrange: %2.1f microns; Yrange: %2.1f',xRange, yRange));


% Plot spectra of eye movements
subplot('Position', subplotPosVectors(figRow,3).v);
hold on;
indices = find(tfAxis < 500);
plot(tfAxis(indices), emSpectrumX(indices), 'm-', 'LineWidth', 1.5);
plot(tfAxis(indices), emSpectrumY(indices), 'b-', 'LineWidth', 1.5);
if (~isempty(emSpectrumXo))
    plot(tfAxis(indices), emSpectrumXo(indices), 'k-', 'LineWidth', 1.5);
    plot(tfAxis(indices), emSpectrumYo(indices), 'k-', 'LineWidth', 1.5);
end
set(gca, 'XLim', [0.1 500], 'YLim', [1e1 1e5], 'XScale', 'log', 'YScale', 'log');
set(gca, 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
if (figRow == 3)
xlabel('frequency (Hz)',  'FontWeight', 'bold');
end
ylabel('power',  'FontWeight', 'bold');
hL = legend({'x-pos', 'y-pos'}, 'Location', 'SouthWest');
set(hL, 'FontSize', 14);
axis square
grid on
title(sprintf('mean spectra of X/Y eye pos movements\n(nTrials: %d)', size(theEMpaths,1)));

drawnow;
end

