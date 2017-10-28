function t_emDifferentConeMosaics
%%t_emDifferentConeMosaics  Tests whether the eye movement spatial/temporal dynamics are independent of the spatial characteristics of the cone mosaic.  
%
% Description:
%    Demonstrates eye movements dynamics for three types of cone mosaics: 
%     (a) a rectangular cone mosaic,
%     (b) a hexagonal cone mosaic with regular cone spacing and larger (than the default values in isetbio) aperture & spacing 
%     (c) a hexagonal cone mosaic with eccentricity-based cone spacing and smaller than the default value cone aperture
%
%    Each of the three rows in the generated figure corresponds to one of the three cone mosaics studied.
%    The first column  plots the cone mosaic with eye movement from individual trial ssuperimposed in red.
%    The second column plots the x- and y-components of eye movement from individual trials.
%    The third column  plots the spatial distribution of x- and y-eye positions across all computed trials (256 here)
%    The fourth column plots the mean spectra of the x- and y-components across all computed trials (256 here)
%
%    Notes: The eye movements have similar spatial and temporal dynamics across all 3 mosaics. 
%           This demonstrates that the eye movement dynamics are independent of the employed mosaic's 
%           spatial characteristics (pattern size, the cone spacing, and cone aperture).
%

% History:
%
% 10/03/17  npc  (c) isetbio team, 2017.

%% Initialize
ieInit;

%% Set the random seed
rng(1);

%% Set parameters
% The mosaic's field of view
fovDegs = 0.15;

% The cone integration time - also the time base for eye movements
integrationTimeMillisecs = 5;

% Number of eye movement trials to generate
nTrials = 256;

% Generate this many eye movements per trial
eyeMovementsPerTrial = 1100;

% Spatial position binning (in microns)
posBinsMicrons = [-28:2:28];

%% Run the default rect mosaic
eccBasedConeDensity = [];
customLamda = []; 
customInnerSegmentDiameter = [];
resamplingFactor = [];
[coneMosaicOBJ0, theEMpathsMicrons0, emSpectrumX0, emSpectrumY0, tfAxis, xPosMean0, yPosMean0, xPosStd0, yPosStd0] = ...
    t_emMosaic('rect', fovDegs, integrationTimeMillisecs, nTrials, eyeMovementsPerTrial, posBinsMicrons, ...
    eccBasedConeDensity, customLamda, customInnerSegmentDiameter, resamplingFactor);

%% Run a regular hex cone mosaic with large separation (customLambda: 7 microns) and large inner segment diameter (4 microns)
eccBasedConeDensity = false; 
customLamda = 7; 
customInnerSegmentDiameter = 4;
resamplingFactor = 6;
[coneMosaicOBJ1, theEMpathsMicrons1, emSpectrumX1, emSpectrumY1, tfAxis, xPosMean1, yPosMean1, xPosStd1, yPosStd1] = ...
    t_emMosaic('hex', fovDegs, integrationTimeMillisecs, nTrials, eyeMovementsPerTrial, posBinsMicrons, ...
    eccBasedConeDensity, customLamda, customInnerSegmentDiameter, resamplingFactor);

%% Run an ecc-based hex cone mosaic with small inner segment diameter (1 micron)
eccBasedConeDensity = true; 
customInnerSegmentDiameter = 1; 
customLamda = [];
[coneMosaicOBJ2, theEMpathsMicrons2, emSpectrumX2, emSpectrumY2, tfAxis, xPosMean2, yPosMean2, xPosStd2, yPosStd2] = ...
    t_emMosaic('hex', fovDegs, integrationTimeMillisecs, nTrials, eyeMovementsPerTrial, posBinsMicrons, ...
    eccBasedConeDensity, customLamda, customInnerSegmentDiameter, resamplingFactor);

%% Display results
for trialToDisplay = 1:10
plotEMs(1, coneMosaicOBJ0, theEMpathsMicrons0, trialToDisplay, emSpectrumX0, emSpectrumY0, [], [], tfAxis, xPosMean0, yPosMean0, xPosStd0, yPosStd0, posBinsMicrons);
plotEMs(2, coneMosaicOBJ1, theEMpathsMicrons1, trialToDisplay, emSpectrumX1, emSpectrumY1, emSpectrumX0, emSpectrumY0, tfAxis, xPosMean1, yPosMean1, xPosStd1, yPosStd1, posBinsMicrons);
plotEMs(3, coneMosaicOBJ2, theEMpathsMicrons2, trialToDisplay, emSpectrumX2, emSpectrumY2, emSpectrumX0, emSpectrumY0, tfAxis, xPosMean2, yPosMean2, xPosStd2, yPosStd2, posBinsMicrons);
drawnow
end

end

% Method to run the rect mosaic simulation
function [cm, theEMpathsMicrons, emSpectrumX, emSpectrumY, tfAxis, ...
    xDistributionMean, yDistributionMean, xDistributionStd, yDistributionStd] = ...
    t_emMosaic(mosaicType, fovDegs, integrationTimeMillisecs, ... ...
    nTrials, eyeMovementsPerTrial, posBinsMicrons, ...
    eccBasedConeDensity, customLamda, customInnerSegmentDiameter, resamplingFactor)

%% Generate mosaic
if (strcmp(mosaicType, 'rect'))
    % Generate a rect mosaic
    cm = coneMosaic();
    cm.setSizeToFOV(fovDegs);            
else
    % Generate a hex mosaic
    cm = coneMosaicHex(resamplingFactor, ...
        'fovDegs', fovDegs, ...
        'eccBasedConeDensity', eccBasedConeDensity, ... 
        'customLambda', customLamda, ...
        'customInnerSegmentDiameter', customInnerSegmentDiameter, ...
        'latticeAdjustmentPositionalToleranceF', 0.1, ...
        'latticeAdjustmentDelaunayToleranceF', 0.1 ...
    );
end

%% Set the integration time - also the time base for eye movements
cm.integrationTime = integrationTimeMillisecs/1000;  

%% Generate eye movement paths
theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
for iTrial= 1:nTrials
    [theEMpaths(iTrial, :,:), theEMpathsMicrons(iTrial, :,:)] = cm.emGenSequence(eyeMovementsPerTrial);
end
theEMpaths = [];

%% Compute spatial stats
[xDistributionMean, yDistributionMean, xDistributionStd, yDistributionStd] = computeEMstats(theEMpathsMicrons, posBinsMicrons);

%% Compute temporal spectra
[emSpectrumX, emSpectrumY, tfAxis] = computeEMspectrum(cm.timeAxis*1000, theEMpathsMicrons);
end

% ----- Spatial stats analyis routine ------
function [xDistributionMean, yDistributionMean, xDistributionStd, yDistributionStd] = computeEMstats(theEMpathsMicrons, posBinsMicrons)
    
    dt = posBinsMicrons(2)-posBinsMicrons(1);
    posBinsMicrons(end+1) = posBinsMicrons(end)+dt;

    for emTrial = 1:size(theEMpathsMicrons,1)
        xPath = squeeze(theEMpathsMicrons(emTrial,:,1));
        yPath = squeeze(theEMpathsMicrons(emTrial,:,2));
        [xDistribution(emTrial,:), xBins] = histcounts(xPath, posBinsMicrons);
        [yDistribution(emTrial,:), yBins] = histcounts(yPath, posBinsMicrons);
    end

    xDistributionMean = mean(xDistribution,1);
    yDistributionMean = mean(yDistribution,1);
    
    xDistributionStd  = std(xDistribution,0,1);
    yDistributionStd  = std(yDistribution,0,1);
end


% ----- Spectrum analysis routine ----
function [emSpectrumX, emSpectrumY, tfAxis] = computeEMspectrum(timeAxisMillisecs, theEMpathsMicrons)
fN = 2^(ceil(log(size(theEMpathsMicrons,2))/log(2)));
for emTrial = 1:size(theEMpathsMicrons,1)
    emSpectrumX(emTrial,:) = abs(fft(squeeze(theEMpathsMicrons(emTrial,:,1)), fN));
    emSpectrumY(emTrial,:) = abs(fft(squeeze(theEMpathsMicrons(emTrial,:,2)), fN));
end
emSpectrumX = fftshift(squeeze(mean(emSpectrumX, 1)));
emSpectrumX = emSpectrumX(fN/2:end);
emSpectrumY = fftshift(squeeze(mean(emSpectrumY, 1)));
emSpectrumY = emSpectrumY(fN/2:end);

maxF = 1000/(2*(timeAxisMillisecs(2)-timeAxisMillisecs(1)));
deltaF = maxF / (fN/2);
tfAxis = (0:fN/2)*deltaF;
end


%----- PLOTTING ROUTINE -------
function plotEMs(figRow, cm, theEMpathsMicrons, iTrial, emSpectrumX, emSpectrumY, emSpectrumXo, emSpectrumYo, tfAxis, xPosMean, yPosMean, xPosStd, yPosStd, posBinsMicrons)

% convert cone positions to microns
if (isa(cm, 'coneMosaicHex'))
    conePositionsMicrons = cm.coneLocs*1e6/cm.resamplingFactor;
else
    conePositionsMicrons = cm.coneLocs*1e6;
end

% convert aperture size to microns
apertureSizeMicrons = cm.pigment.width*1e6;
spatialExtent = [posBinsMicrons(1) posBinsMicrons(end)];
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
    set(hFig, 'Position', [10 10 1700 1300], 'Color', [1 1 1]);
end

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', 4, ...
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
        'overlayEMpathMicrons', squeeze(theEMpathsMicrons(iTrial,:,:)), ...  % expects emPath in microns, not cone units
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
    title(sprintf('Hex mosaic,  pattern sample size: %2.1f um\ncone sep: %g um, aperture: %2.1f um (trial #%d)', cm.patternSampleSize(1)*1e6, customLambda, cm.pigment.pdWidth*1e6, iTrial));
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
    title(sprintf('Rect mosaic, pattern sample size: %2.1f um\ncone sep: %g um, aperture: %2.1f um (trial #%d)', cm.patternSampleSize(1)*1e6, cm.patternSampleSize(1)*1e6, cm.pigment.pdWidth*1e6, iTrial));
end
ylabel('space (microns)',  'FontWeight', 'bold');
if (figRow == 3)
    xlabel('space (microns)',  'FontWeight', 'bold');
end
set(gca, 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
box(gca, 'on'); grid(gca, 'off');
    

% Plot the time course of the x and y eye movement
subplot('Position', subplotPosVectors(figRow,2).v);
plot(timeAxisMillisecs, theEMpathsMicrons(iTrial,:,1), '-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerEdgeColor', 'm', 'Color', 'm');
hold on;
plot(timeAxisMillisecs, theEMpathsMicrons(iTrial,:,2), '-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'Color', 'b');
set(gca, 'XLim', [timeAxisMillisecs(1) timeAxisMillisecs(end)], 'YLim', spatialExtent, 'YTick', yTicks);
grid on; box off;
set(gca, 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
hL = legend({'x-pos', 'y-pos'}, 'Location', 'SouthWest');
set(hL, 'FontSize', 14);
ylabel('space (microns)',  'FontWeight', 'bold'); 
if (figRow == 3)
    xlabel('time (msec)', 'FontWeight', 'bold'); 
end
title(sprintf('X/Y eye position sequence\n(trial #%d)', iTrial));


subplot('Position', subplotPosVectors(figRow,3).v);
maxCount = max([max(xPosMean+xPosStd) max(yPosMean+yPosStd)]);
hold on
stairs(posBinsMicrons, xPosMean, 'LineWidth', 1.5, 'Color', 'm');
stairs(posBinsMicrons, -yPosMean, 'LineWidth', 1.5, 'Color', 'b');
stairs(posBinsMicrons, xPosMean+xPosStd, 'LineWidth', 1.0, 'LineStyle', ':', 'Color', 'm');
stairs(posBinsMicrons, xPosMean-xPosStd, 'LineWidth', 1.0, 'LineStyle', ':', 'Color', 'm');
stairs(posBinsMicrons, -yPosMean+yPosStd, 'LineWidth', 1.0, 'LineStyle', ':', 'Color', 'b');
stairs(posBinsMicrons, -yPosMean-yPosStd, 'LineWidth', 1.0, 'LineStyle', ':', 'Color', 'b');
plot([posBinsMicrons(1) posBinsMicrons(end)], [0 0], 'k-');
set(gca, 'XLim', [posBinsMicrons(1) posBinsMicrons(end)], 'YLim', maxCount*[-1 1], 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
grid on; box on;
ylabel('power',  'FontWeight', 'bold');
hL = legend({'x-pos', 'y-pos'}, 'Location', 'SouthWest');
set(hL, 'FontSize', 14);
if (figRow == 3)
xlabel('position (um)', 'FontWeight', 'bold');
end
ylabel('counts', 'FontWeight', 'bold');
axis 'square';
title(sprintf('mean and std of X/Y eye positions \n(nTrials: %d)', size(theEMpathsMicrons,1)));
    



% Plot spectra of eye movements
subplot('Position', subplotPosVectors(figRow,4).v);
hold on;
plot(tfAxis, emSpectrumX, 'm-', 'LineWidth', 1.5);
plot(tfAxis, emSpectrumY, 'b-', 'LineWidth', 1.5);
hL = legend({'x-pos', 'y-pos'}, 'Location', 'SouthWest');
if (~isempty(emSpectrumXo))
    plot(tfAxis, emSpectrumXo, 'k--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
    plot(tfAxis, emSpectrumYo, 'k--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.0);
    hL = legend({'x-pos', 'y-pos', 'x-pos (rect mosaic)', 'y-pos (rect mosaic)'}, 'Location', 'SouthWest');
end
set(gca, 'XLim', [0.1 tfAxis(end)], 'YLim', [1e1 1e5], 'XScale', 'log', 'YScale', 'log');
set(gca, 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
if (figRow == 3)
xlabel('frequency (Hz)',  'FontWeight', 'bold');
end
ylabel('power',  'FontWeight', 'bold');
set(hL, 'FontSize', 14, 'Color', [0.9 0.9 0.9]);
axis square
grid on
title(sprintf('mean spectra of X/Y eye pos movements\n(nTrials: %d)', size(theEMpathsMicrons,1)));
end



