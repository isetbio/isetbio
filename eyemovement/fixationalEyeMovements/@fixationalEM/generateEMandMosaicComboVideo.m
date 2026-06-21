function generateEMandMosaicComboVideo(fixEMobj, coneMosaic, varargin)
% Generate a video showing the eye movement on top of the cone mosaic.
%
% Syntax:
%   generateEMandMosaicComboVideo(fixEMobj, coneMosaic, [varargin])
%   fixEMobj.generateEMandMosaicComboVideo(coneMosaic, [varargin])
%
% Description:
%    Generate a video showing the eye movement on top of the cone mosaic.
%
% Inputs:
%    fixEMobj            - Object. A fixationalEM object.
%    coneMosaic          - Object. A cMosaic object.
%    varargin            - (Optional) Any additional arguments required,
%                          see key/values section for more information.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    'visualizedFOVdegs' - Numeric. Visualized FoV in degrees. Default 0.4
%    'displayCrosshairs' - Boolean. Should display Crosshairs? Default true
%    'showMovingMosaicOnSeparateSubFig'
%                        - Boolean. Should we display the mosaic on a
%                          separate sub-figure? Default true.
%    'videoFileName'     - String. The filename to save generated video as.
%                          Default emMosaicCombo.mp4
%    'frameStep'         - Positive integer. Capture every frameStep-th eye
%                          movement sample. The final sample is always
%                          captured. Default 1.
%    'frameRate'         - Positive scalar. Output video frame rate in
%                          frames/second. Default 60.
%    'videoQuality'      - Scalar in [0,100]. MPEG-4 encoding quality.
%                          Default 90.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments
%    06/21/26  BW   Reuse graphics objects, add frame subsampling, and
%                   ensure pixel-safe figure sizing and resource cleanup.
%

p = inputParser;
p.addRequired('fixEMobj', @(x)(isa(x, 'fixationalEM')));
p.addRequired('coneMosaic', @(x)(isa(x, 'cMosaic')));
p.addParameter('visualizedFOVdegs', 1.0, ...
    @(x)(isnumeric(x) && isscalar(x) && isfinite(x) && (x > 0)));
p.addParameter('displayCrosshairs', true, @islogical);
p.addParameter('showMovingMosaicOnSeparateSubFig', true, @islogical);
p.addParameter('videoFileName', 'emMosaicCombo.mp4', ...
    @(x)(ischar(x) || (isstring(x) && isscalar(x))));
p.addParameter('frameStep', 1, ...
    @(x)(isnumeric(x) && isscalar(x) && isfinite(x) && ...
    (x >= 1) && (round(x) == x)));
p.addParameter('frameRate', 60, ...
    @(x)(isnumeric(x) && isscalar(x) && isfinite(x) && (x > 0)));
p.addParameter('videoQuality', 90, ...
    @(x)(isnumeric(x) && isscalar(x) && isfinite(x) && ...
    (x >= 0) && (x <= 100)));
p.parse(fixEMobj, coneMosaic, varargin{:});

visualizedFOVdegs = p.Results.visualizedFOVdegs;
videoFileName = p.Results.videoFileName;
showMovingMosaicOnSeparateSubFig = ...
    p.Results.showMovingMosaicOnSeparateSubFig;
displayCrossHairs = p.Results.displayCrosshairs;
frameStep = p.Results.frameStep;
frameRate = p.Results.frameRate;
videoQuality = p.Results.videoQuality;

if isempty(fixEMobj.emPosMicrons)
    fprintf('No eye movements found in the passed fixationalEM object\n');
    return;
end

nTrials = size(fixEMobj.emPosMicrons, 1);
crossLengthDegs = visualizedFOVdegs / 20;


ticks = linspace(-0.5*visualizedFOVdegs, 0.5*visualizedFOVdegs, 5);
tickLabels = compose('%0.2f', ticks);

% Subfig layout
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 1 + showMovingMosaicOnSeparateSubFig, ...
       'heightMargin', 0.08, ...
       'widthMargin', 0.04, ...
       'leftMargin', 0.06, ...
       'rightMargin', 0.02, ...
       'bottomMargin', 0.08, ...
       'topMargin', 0.02);

if showMovingMosaicOnSeparateSubFig
    figSize = [1100 580];
else
    figSize = [510 580];
end

% Create one fixed-size figure for the entire video. Explicit pixel units
% prevent environmental figure defaults from interpreting figSize as a
% normalized position. Resize is disabled so all video frames have exactly
% the same dimensions.
hFig = figure();
originalFigureUnits = hFig.Units;
hFig.Units = 'pixels';
set(hFig, ...
    'Position', [10 10 figSize(1) figSize(2)], ...
    'Color', [1 1 1], ...
    'Resize', 'off');
hFig.Units = originalFigureUnits;

% Create the axes once. Repeatedly creating figures, axes, and line objects
% puts unnecessary pressure on MATLAB's graphics scene.
if showMovingMosaicOnSeparateSubFig
    axMosaic = axes('Parent', hFig, ...
        'Position', subplotPosVectors(1, 1).v);
    axEMpath = axes('Parent', hFig, ...
        'Position', subplotPosVectors(1, 2).v);
else
    axMosaic = gobjects(0);
    axEMpath = axes('Parent', hFig, ...
        'Position', subplotPosVectors(1, 1).v);
end

% Render the cone mosaic once. Moving the viewport requires only axes-limit
% updates during the frame loop.
if showMovingMosaicOnSeparateSubFig
    coneMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', axMosaic);
    set(axMosaic, 'LineWidth', 1.0, 'FontSize', 20);
    box(axMosaic, 'on');
    mosaicTitle = title(axMosaic, '');
else
    mosaicTitle = gobjects(0);
end

% Configure the eye-movement axes and create persistent graphics handles.
axis(axEMpath, 'xy');
axis(axEMpath, 'image');
set(axEMpath, ...
    'XLim', visualizedFOVdegs * [-0.5 0.5] * 1.05, ...
    'YLim', visualizedFOVdegs * [-0.5 0.5] * 1.05, ...
    'XTick', ticks, ...
    'YTick', ticks, ...
    'XTickLabel', tickLabels, ...
    'YTickLabel', {}, ...
    'LineWidth', 1.0, ...
    'FontSize', 20);
xlabel(axEMpath, 'space (degs)');
emPathTitle = title(axEMpath, '');
xtickangle(axEMpath, 0);
grid(axEMpath, 'on');
box(axEMpath, 'on');
hold(axEMpath, 'on');

axisExtent = visualizedFOVdegs * [-0.5 0.5] * 1.05;
plot(axEMpath, [0 0], axisExtent, '-', 'Color', [0.3 0.3 0.3]);
plot(axEMpath, axisExtent, [0 0], '-', 'Color', [0.3 0.3 0.3]);

previousTrialLines = gobjects(nTrials, 1);
for k = 1:nTrials
    previousTrialLines(k) = plot(axEMpath, nan, nan, 'k-', ...
        'LineWidth', 2.0);
end
currentPathHalo = plot(axEMpath, nan, nan, '-', ...
    'Color', [1 0.5 0.5], 'LineWidth', 3.0);
currentPathLine = plot(axEMpath, nan, nan, 'r-', 'LineWidth', 1.5);
verticalCrosshair = plot(axEMpath, nan, nan, 'k-', 'LineWidth', 1.5);
horizontalCrosshair = plot(axEMpath, nan, nan, 'k-', 'LineWidth', 1.5);
hold(axEMpath, 'off');

if displayCrossHairs
    crosshairVisibility = 'on';
else
    crosshairVisibility = 'off';
end
set([verticalCrosshair horizontalCrosshair], ...
    'Visible', crosshairVisibility);

% Open the video stream after graphics initialization. The cleanup object
% closes both the stream and figure if rendering or encoding throws an error.
videoOBJ = VideoWriter(videoFileName, 'MPEG-4');
videoOBJ.FrameRate = frameRate;
videoOBJ.Quality = videoQuality;
videoOBJ.open();
resourceCleanup = onCleanup(@() cleanupVideoResources(videoOBJ, hFig));

for iTrial = 1:nTrials
    theEMpathDegs = squeeze(fixEMobj.emPosArcMin(iTrial, :, :))/60;

    % Display completed paths from earlier trials. Later trials remain hidden
    % until they become part of the accumulated history.
    for k = 1:nTrials
        if k < iTrial
            previousPathDegs = squeeze(fixEMobj.emPosArcMin(k, :, :))/60;
            set(previousTrialLines(k), ...
                'XData', previousPathDegs(:, 1), ...
                'YData', previousPathDegs(:, 2));
        else
            set(previousTrialLines(k), 'XData', nan, 'YData', nan);
        end
    end

    set([currentPathHalo currentPathLine verticalCrosshair horizontalCrosshair], ...
        'XData', nan, 'YData', nan);

    % Capture every frameStep-th sample and always include the final sample.
    nTimeBins = size(theEMpathDegs, 1);
    capturedTimeBins = unique([1:frameStep:nTimeBins nTimeBins]);

    for timeBin = capturedTimeBins
        currentX = theEMpathDegs(1:timeBin, 1);
        currentY = theEMpathDegs(1:timeBin, 2);
        set([currentPathHalo currentPathLine], ...
            'XData', currentX, 'YData', currentY);

        xo = theEMpathDegs(timeBin, 1);
        yo = theEMpathDegs(timeBin, 2);

        if displayCrossHairs
            set(verticalCrosshair, ...
                'XData', xo + [0 0], ...
                'YData', yo + crossLengthDegs * [-1 1]);
            set(horizontalCrosshair, ...
                'XData', xo + crossLengthDegs * [-1 1], ...
                'YData', yo + [0 0]);
        end

        if showMovingMosaicOnSeparateSubFig
            set(axMosaic, ...
                'XLim', coneMosaic.eccentricityDegs(1) + xo + ...
                    visualizedFOVdegs * [-0.5 0.5] * 1.05, ...
                'YLim', coneMosaic.eccentricityDegs(2) + yo + ...
                    visualizedFOVdegs * [-0.5 0.5] * 1.05, ...
                'XTick', coneMosaic.eccentricityDegs(1) + xo + ...
                    visualizedFOVdegs * [-0.5 0.0 0.5], ...
                'YTick', coneMosaic.eccentricityDegs(2) + yo + ...
                    visualizedFOVdegs * [-0.5 0.0 0.5]);
            set(mosaicTitle, 'String', ...
                sprintf('%2.1f msec', 1000 * fixEMobj.timeAxis(timeBin)));
        end

        if nTrials > 1
            if showMovingMosaicOnSeparateSubFig
                titleString = sprintf('trial #%d', iTrial);
            else
                titleString = sprintf('%2.1f msec (trial #%d)', ...
                    1000 * fixEMobj.timeAxis(timeBin), iTrial);
            end
        else
            titleString = sprintf('%2.1f msec', ...
                1000 * fixEMobj.timeAxis(timeBin));
        end
        set(emPathTitle, 'String', titleString);

        % getframe performs the required full graphics update, so a separate
        % drawnow call would duplicate work.
        videoOBJ.writeVideo(getframe(hFig));
    end
end

videoOBJ.close();
fprintf('File saved in %s\n', videoFileName);

% Invoke cleanup now so the large mosaic figure does not remain resident.
clear resourceCleanup;

end

function cleanupVideoResources(videoOBJ, hFig)
% Close resources during normal completion and after MATLAB-level errors.

try
    videoOBJ.close();
catch
    % The writer may not have opened or may already be closed.
end

if isgraphics(hFig)
    close(hFig);
end
end
