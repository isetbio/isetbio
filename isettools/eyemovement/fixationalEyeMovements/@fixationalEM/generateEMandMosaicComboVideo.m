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
%    coneMosaic          - Object. A coneMosaicHex object.
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
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments
%

p = inputParser;
p.addRequired('fixEMobj', @(x)(isa(x, 'fixationalEM')));
p.addRequired('coneMosaic', @(x)(isa(x, 'coneMosaicHex')));
p.addParameter('visualizedFOVdegs', 0.4, @isnumeric);
p.addParameter('displayCrosshairs', true, @islogical);
p.addParameter('showMovingMosaicOnSeparateSubFig', true, @islogical);
p.addParameter('videoFileName', 'emMosaicCombo.mp4', @ischar);
p.parse(fixEMobj, coneMosaic, varargin{:});

visualizedFOVdegs = p.Results.visualizedFOVdegs;
videoFileName = p.Results.videoFileName;
showMovingMosaicOnSeparateSubFig = ...
    p.Results.showMovingMosaicOnSeparateSubFig;
displayCrossHairs = p.Results.displayCrosshairs;

if isempty(fixEMobj.emPosMicrons)
    fprintf('No eye movements found in the passed fixationalEM object\n');
    return;
end

nTrials = size(fixEMobj.emPosMicrons, 1);
crossLengthDegs = visualizedFOVdegs / 20;
crossLengthMeters = crossLengthDegs * coneMosaic.micronsPerDegree * 1e-6;

visualizedSpaceMeters = ...
    visualizedFOVdegs * coneMosaic.micronsPerDegree * 1e-6;
ticksDegs = -0.5:0.1:0.5;
ticks = ticksDegs * coneMosaic.micronsPerDegree * 1e-6;
tickLabels = sprintf('%2.1f\n', ticksDegs);

% Subfig layout
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 1, ...
       'colsNum', 1 + showMovingMosaicOnSeparateSubFig, ...
       'heightMargin', 0.08, ...
       'widthMargin', 0.04, ...
       'leftMargin', 0.03, ...
       'rightMargin', 0.02, ...
       'bottomMargin', 0.06, ...
       'topMargin', 0.02);

if (showMovingMosaicOnSeparateSubFig)
    figSize = [1100 580];
else
    figSize = [510 580];
end

% Open video stream
videoOBJ = VideoWriter(videoFileName, 'MPEG-4'); % H264 format
videoOBJ.FrameRate = 60;
videoOBJ.Quality = 100;
videoOBJ.open();

for iTrial = 1:nTrials
    % Make figure
    hFig = figure();
    
    set(hFig, 'Position', [10 10 figSize(1) figSize(2)], 'Color', [1 1 1]);
    % Make subplots axes
    if (showMovingMosaicOnSeparateSubFig)
        axMosaic = subplot('Position', subplotPosVectors(1, 1).v);
        axEMpath = subplot('Position', subplotPosVectors(1, 2).v);
    else
        axEMpath = subplot('Position', subplotPosVectors(1, 1).v);
    end

    % Render the mosaic
    if (showMovingMosaicOnSeparateSubFig)
        coneMosaic.visualizeGrid('axeshandle', axMosaic, ...
            'backgroundColor', [.75 .75 .75], ...
            'foregroundColor', [0 0 0], ...
            'visualizedConeAperture', 'geometricArea', ...
            'labelConeTypes', true, ...
            'apertureshape', 'disks');
        hold(axMosaic, 'on');
        if (displayCrossHairs)
            plot(axMosaic, [0 0], crossLengthMeters * [-1 1], 'k-', ...
                'LineWidth', 2.0);
            plot(axMosaic, crossLengthMeters * [-1 1], [0 0], 'k-', ...
                'LineWidth', 2.0);
        end
        set(axMosaic, 'LineWidth', 1.0, 'FontSize', 20);
        hold(axMosaic, 'off');
        box(axMosaic, 'on');
    end

    
    
    % emPath for this trial
    theEMpathMeters = squeeze(fixEMobj.emPosMicrons(iTrial, :, :)) * 1e-6;

    % Render frames of the combo video, one frame per each timeBin
    for timeBin = 1:size(theEMpathMeters, 1)
        % Move the mosaic
        if (showMovingMosaicOnSeparateSubFig)
            xo = theEMpathMeters(timeBin, 1);
            yo = theEMpathMeters(timeBin, 2);
            set(axMosaic, 'XLim', xo + visualizedSpaceMeters * ...
                [-0.5 0.5] * 1.05, 'YLim', yo + visualizedSpaceMeters * ...
                [-0.5 0.5] * 1.05, 'XTickLabel', {}, 'YTickLabel', {});
            set(axMosaic, 'FontSize', 20, 'LineWidth', 1.0);
            title(axMosaic, ...
                sprintf('%2.1f msec', 1000 * fixEMobj.timeAxis(timeBin)));
        end

        
        % plot emPaths for previous trials
        for k = 1: iTrial-1
            theEMpathMetersOld = squeeze(fixEMobj.emPosMicrons(k, :, :)) * 1e-6;
            plot(axEMpath, theEMpathMetersOld(:, 1), ...
                theEMpathMetersOld(:, 2), 'k-',  ...
                'LineWidth', 2.0);
            if (k == 1)
                hold(axEMpath, 'on');
            end
        end
        
        
        
        % Update the emPath
        plot(axEMpath, theEMpathMeters(1:timeBin, 1), ...
            theEMpathMeters(1:timeBin, 2), 'r-', 'Color', [1 0.5 0.5], ...
            'LineWidth', 3.0);
        hold(axEMpath, 'on');
        plot(axEMpath, theEMpathMeters(1:timeBin, 1), ...
            theEMpathMeters(1:timeBin, 2), 'r-', 'LineWidth', 1.5);

        if (displayCrossHairs)
            plot(axEMpath, theEMpathMeters(timeBin, 1) + [0 0], ...
                theEMpathMeters(timeBin, 2) + crossLengthMeters * ...
                [-1 1], 'k-', 'LineWidth', 1.5);
            plot(axEMpath, theEMpathMeters(timeBin, 1) + ...
                crossLengthMeters * [-1 1], theEMpathMeters(timeBin, 2) ...
                + [0 0], 'k-', 'LineWidth', 1.5);
        end

        axis(axEMpath, 'xy');
        axis(axEMpath, 'image');
        plot(axEMpath, [0 0], 1e-6 * [-300 300], '-', ...
            'Color', [0.3 0.3 0.3]);
        plot(axEMpath, 1e-6 * [-300 300], [0 0], '-', ...
            'Color', [0.3 0.3 0.3]);
        hold(axEMpath, 'off');

        % Set ticks and limits
        set(axEMpath, 'XLim', visualizedSpaceMeters * [-0.5 0.5] * ...
            1.05, 'YLim', visualizedSpaceMeters * [-0.5 0.5] * 1.05, ...
            'XTick', ticks, 'YTick', ticks, ...
            'XTickLabels', tickLabels, 'YTickLabels', {}, 'FontSize', 16);
        xlabel(axEMpath, 'space (degs)');
        grid(axEMpath, 'on');
        box(axEMpath, 'on');
        set(axEMpath, 'LineWidth', 1.0, 'FontSize', 20);

        if (nTrials > 1)
            if (~showMovingMosaicOnSeparateSubFig)
                title(axEMpath, sprintf('%2.1f msec (trial #%d)', ...
                    1000 * fixEMobj.timeAxis(timeBin), iTrial));
            else
                title(axEMpath, sprintf('trial #%d', iTrial));
            end
        else
            title(axEMpath, sprintf('%2.1f msec', ...
                1000 * fixEMobj.timeAxis(timeBin)));
        end

        % Render
        drawnow;

        % Add video frame
        videoOBJ.writeVideo(getframe(hFig));
    end % for k
end % for iTrial

% Close video stream
videoOBJ.close();
fprintf('File saved in %s\n', videoFileName);

end
