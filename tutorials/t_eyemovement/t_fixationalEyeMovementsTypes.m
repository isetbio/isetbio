function t_fixationalEyeMovementsTypes(varargin)
% Examine eye movement paths using different 'microSaccadeType' parameters.
%
% Syntax:
%   t_fixationalEyeMovementsTypes
%
% Description:
%    The fixation maps are a bit different for microsaccades that are
%    generated using the 'heatmap/fixation based' strategy vs. the 'stats
%    based' strategy. The heatmap/fixation strategy results in fixation
%    maps that are a bit wider along the horizontal and vertical axes.
%
% Inputs:
%    None required.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    'random seed'              - Numeric. Seed for the rng. Default 1
%    'n trials'                 - Numeric. Number of trials. Default 100.
%    'visualized single trials' - Numeric. The number of trials to display
%                                 in addition to the compilation, and the
%                                 fixation map. Default is 3.
%    'sample time seconds'      - Numeric. The duration of the sample in
%                                 seconds. Default is 0.001.
%    'em duration seconds'      - Numeric. The duration of the eye
%                                 movement in seconds. Default is 0.050.
%    'compute velocity'         - Boolean. Whether or not to compute the
%                                 eye movement velocity. Default is false.
%    'use par for'              - Boolean. Whether to compute using a
%                                 parfor loop (to take advantage of
%                                 multiple processors) of to use a single
%                                 processor. Default false, i.e., use a
%                                 single processor.
%

% History
%    02/06/18  npc  Wrote it.
%    02/07/18  npc  Comments.
%    07/09/18  jnm  Formatting. Written as a function, so formatted as one.

% Examples:
%{
    t_fixationalEyeMovementsTypes('random seed', 1);
    t_fixationalEyeMovementsTypes('random seed', 1, 'n trials', 10);
    t_fixationalEyeMovementsTypes('random seed', 1, 'n trials', 10, ...
        'visualizedsingletrials', 2);
%}

%%
p = inputParser;
varargin = ieParamFormat(varargin);

p.addParameter('ntrials', 100, @isscalar);
p.addParameter('sampletimeseconds', 0.001, @isscalar)
p.addParameter('emdurationseconds', 0.050, @isscalar);
p.addParameter('computevelocity', false, @islogical);
p.addParameter('randomseed', 1, @isscalar);
p.addParameter('visualizedsingletrials', 3, @isscalar);
p.addParameter('useparfor', false, @islogical);

p.parse(varargin{:});

emDurationSeconds = p.Results.emdurationseconds;
sampleTimeSeconds = p.Results.sampletimeseconds;
nTrials = p.Results.ntrials;
computeVelocity = p.Results.computevelocity;
randomSeed = p.Results.randomseed;
useParfor = p.Results.useparfor;
visualizedSingleTrials = p.Results.visualizedsingletrials;

%% Initialize object
fixEMobj = fixationalEM();

% First case: No micro-saccades, only drift
% Set all params to their default value
fixEMobj.setDefaultParams();
fixEMobj.microSaccadeType = 'none';
fixEMobj.randomSeed = randomSeed;
fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, ...
    computeVelocity, 'useParfor', useParfor);
plotTrials(fixEMobj, 1, visualizedSingleTrials); drawnow;

% Second case: 'stats based' micro-saccades
fixEMobj.setDefaultParams();
fixEMobj.microSaccadeType = 'stats based';
fixEMobj.randomSeed = randomSeed;
fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, ...
    computeVelocity, 'useParfor', useParfor);
plotTrials(fixEMobj, 2, visualizedSingleTrials); drawnow;

% Third case: 'heatmap/fixation based' micro-saccades
%
% Illustrate use of a custom seed for the rng as well.
fixEMobj.setDefaultParams();
fixEMobj.randomSeed = 678;
fixEMobj.microSaccadeType = 'heatmap/fixation based';
fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, ...
    computeVelocity, 'useParfor', useParfor);
plotTrials(fixEMobj, 3, visualizedSingleTrials);

end

function plotTrials(fixEMobj, rowNo, visualizedSingleTrials)
% Plot the fixational eye movement trials
%
% Syntax:
%   plotTrials(fixEMobj, rowNo, visualizedSingleTrials)
%
% Description:
%    Plot the fixational eye movements' trials using the provided
%    parameters to control the number of trials displayed.
%
% Inputs:
%    fixEMobj               - Object. A fixational eye movement object.
%    rowNo                  - Numeric. The row index.
%    visualizedSingleTrials - Numeric. The number of trials to display
%                             (Final display adds an additional two to this
%                             number for other visualizations.)
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

nTrials = size(fixEMobj.emPosArcMin, 1);

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', 3, ...
    'colsNum', visualizedSingleTrials + 2, ...
    'heightMargin', 0.01, ...
    'widthMargin', 0.01, ...
    'leftMargin', 0.01, ...
    'rightMargin', 0.00, ...
    'bottomMargin', 0.01, ...
    'topMargin', 0.01);

% Plot all trials
xyRange = [-20 20];
tickLabel = xyRange(1):5:xyRange(2);

if (rowNo == 1)
    hFig = figure(1);
    clf;
    set(hFig, 'Position', [10 10 1600 1100]);
else
    figure(1);
end

for iTrial = 1:visualizedSingleTrials + 2
    subplot('Position', subplotPosVectors(rowNo, iTrial).v);
    if (iTrial <= visualizedSingleTrials)
        plot(squeeze(fixEMobj.emPosArcMin(iTrial, :, 1)), ...
            squeeze(fixEMobj.emPosArcMin(iTrial, :, 2)), 'k-');
        title(sprintf('trial: %d', iTrial));
    elseif (iTrial == visualizedSingleTrials + 1)
        hold on
        for k = 1:nTrials
            plot(squeeze(fixEMobj.emPosArcMin(k, :, 1)), ...
                squeeze(fixEMobj.emPosArcMin(k, :, 2)), 'k-');
        end
        title(sprintf('%d trials\nsaccade type:\n''%s''', ...
            nTrials, fixEMobj.microSaccadeType));
    else
        binWidthArcMin = 0.5;
        [fixationMap, fixationMapSupportX, fixationMapSupportY, ...
            fixationMapXSlice, fixationMapYSlice] = ...
            fixEMobj.computeFixationMap(fixEMobj.timeAxis, ...
            fixEMobj.emPosArcMin, xyRange, binWidthArcMin);

        contourf(fixationMapSupportX, fixationMapSupportY, fixationMap, ...
            0:0.05:1, 'LineColor', [.5 0.5 0.5]);
        hold on;
        plot(fixationMapSupportX, ...
            xyRange(1) + fixationMapXSlice * xyRange(2) * 0.9, '-', ...
            'Color', [1 0 0], 'LineWidth', 1.5);
        plot(xyRange(1) + fixationMapYSlice * xyRange(2) * 0.9, ...
            fixationMapSupportY, '-', 'Color', [0 0 1], 'LineWidth', 1.5);
        title('Fixation Map');
    end
    hold on
    plot(xyRange, xyRange * 0, 'k-');
    plot(xyRange * 0, xyRange, 'k-');
    set(gca, 'XLim', xyRange, 'YLim', xyRange, 'XTick', tickLabel, ...
        'YTick', tickLabel, 'XTickLabel', ...
        sprintf('%2.1f\n', tickLabel), 'YTickLabel', {});
    if (iTrial<visualizedSingleTrials + 2)
        set(gca, 'XTickLabel', {});
    else
        xlabel('arc min');
    end
    grid on;
    box on;
    axis 'square';
    axis 'xy'
end

colormap(brewermap(1024, 'Greys'));
end
