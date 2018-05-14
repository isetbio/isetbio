function t_fixationalEyeMovementsTypes(varargin)
% Examine eye movement paths produced by different values of the
% 'microSaccadeType' parameter.
% 
% Note that the fixation maps are a bit different for microsaccades that
% are generated using the 'heatmap/fixation based' strategy vs. the
% 'stats based' strategy, with the former strategy resulting in fixation
% maps that are a bit wider along the horizontal and vertical axes.
% 

% History
%   02/06/18  npc  Wrote it.
%   02/07/18  npc  Comments.

    p = inputParser;
    varargin = ieParamFormat(varargin);
    
    p.addParameter('ntrials',100,@isscalar);
    p.addParameter('sampletimeseconds',0.001,@isscalar)
    p.addParameter('emdurationseconds',0.050,@isscalar);
    p.addParameter('useparfor',true,@islogical);

    p.parse(varargin{:});
    
    emDurationSeconds = p.Results.emdurationseconds;
    sampleTimeSeconds = p.Results.sampletimeseconds;
    nTrials           = p.Results.ntrials;
    useParfor         = p.Results.useparfor;
    
    % Generate eye movement data for 2048 1-second long trials, with a
    % sample time of 1 millisecond.
    % emDurationSeconds = 1.0; sampleTimeSeconds = 1/1000; nTrials = 2048;
    
    % Do not compute velocity of eye movements
    computeVelocity = false;
    
    % Random seed to be used in all eye movement compute() calls
    randomSeed = 1;
    
    % Visualize the 5 first trials
    visualizedSingleTrials = 3;
    
    % Initialize object
    fixEMobj = fixationalEM();
    
    % First case: No micro-saccades, only drift
    % Set all params to their default value
    fixEMobj.setDefaultParams();
    fixEMobj.microSaccadeType = 'none';
    fixEMobj.randomSeed = randomSeed;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, ...
        computeVelocity, 'useParfor', useParfor);
    plotTrials(fixEMobj, 1, visualizedSingleTrials);
    
    % Second case: 'stats based' micro-saccades
    fixEMobj.setDefaultParams();
    fixEMobj.microSaccadeType = 'stats based';
    fixEMobj.randomSeed = randomSeed;
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, ...
        computeVelocity, 'useParfor', useParfor);
    plotTrials(fixEMobj, 2, visualizedSingleTrials);
    
    % Third case: 'heatmap/fixation based' micro-saccades
    fixEMobj.setDefaultParams();
    fixEMobj.randomSeed = 678;
    fixEMobj.microSaccadeType = 'heatmap/fixation based';
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, ...
        computeVelocity, 'useParfor', useParfor);
    plotTrials(fixEMobj, 3, visualizedSingleTrials);

end

function plotTrials(fixEMobj, rowNo, visualizedSingleTrials)
    
    nTrials = size(fixEMobj.emPosArcMin,1);
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', visualizedSingleTrials+2, ...
           'heightMargin',   0.01, ...
           'widthMargin',    0.01, ...
           'leftMargin',     0.01, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.01, ...
           'topMargin',      0.01);
       
    % Plot all trials
    xyRange = [-20 20];
    tickLabel = [xyRange(1):5:xyRange(2)];
     
    if (rowNo == 1)
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1600 900]);
    else
        figure(1);
    end
    
    for iTrial = 1:visualizedSingleTrials+2
        subplot('Position', subplotPosVectors(rowNo,iTrial).v);
        if (iTrial <= visualizedSingleTrials)
            plot(squeeze(fixEMobj.emPosArcMin(iTrial,:,1)), squeeze(fixEMobj.emPosArcMin(iTrial,:,2)), 'k-');
            title(sprintf('trial: %d', iTrial));
        elseif (iTrial == visualizedSingleTrials+1)
            hold on
            for k = 1:nTrials
                plot(squeeze(fixEMobj.emPosArcMin(k,:,1)), squeeze(fixEMobj.emPosArcMin(k,:,2)), 'k-');
            end
            title(sprintf('%d trials\nsaccade type:\n''%s''', nTrials,fixEMobj.microSaccadeType));
        else
            binWidthArcMin = 0.5;
            [fixationMap, fixationMapSupportX, fixationMapSupportY, fixationMapXSlice, fixationMapYSlice] = ...
                fixEMobj.computeFixationMap(fixEMobj.timeAxis, fixEMobj.emPosArcMin, ...
                    xyRange, binWidthArcMin);
            
            contourf(fixationMapSupportX, fixationMapSupportY, fixationMap, 0:0.05:1, 'LineColor', [.5 0.5 0.5]); hold on; 
            plot(fixationMapSupportX, xyRange(1)+fixationMapXSlice*xyRange(2)*0.9, '-', 'Color', [1 0 0], 'LineWidth', 1.5);
            plot(xyRange(1)+fixationMapYSlice*xyRange(2)*0.9, fixationMapSupportY, '-', 'Color', [0 0 1], 'LineWidth', 1.5);
            title('Fixation Map');
        end
        hold on
        plot(xyRange, xyRange*0, 'k-');
        plot(xyRange*0, xyRange, 'k-');
        set(gca, 'XLim', xyRange, 'YLim', xyRange, 'XTick', tickLabel, 'YTick', tickLabel, ...
                'XTickLabel', sprintf('%2.1f\n', tickLabel), 'YTickLabel', {});
        if (iTrial<visualizedSingleTrials+2)
            set(gca, 'XTickLabel', {});
        else
            xlabel('arc min');
        end
        grid on; box on; axis 'square'; axis 'xy'
    end
    colormap(brewermap(1024, 'Greys'));
end
