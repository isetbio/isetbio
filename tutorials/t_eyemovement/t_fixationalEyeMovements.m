function t_fixationalEyeMovements
% Examine eye movement paths produced by different values of the
% 'microSaccadeType' parameter.
%

% History
%   02/06/18  npc  Wrote it.

    close all;
    
    emDurationSeconds = 1.0; sampleTimeSeconds = 1/1000; nTrials = 20;
    
    % Initialize object
    fixEMobj = fixationalEM();
    computeVelocity = false;
    
    fixEMobj.setDefaultParams();
    fixEMobj.randomSeed = 678;
    fixEMobj.microSaccadeType = 'none';
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, computeVelocity);
    plotAllTrials(fixEMobj, 1);
    
    fixEMobj.setDefaultParams();
    fixEMobj.randomSeed = 678;
    fixEMobj.microSaccadeType = 'stats based';
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, computeVelocity);
    plotAllTrials(fixEMobj, 2);
    
    fixEMobj.setDefaultParams();
    fixEMobj.randomSeed = 678;
    fixEMobj.microSaccadeType = 'heatmap/fixation based';
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials, computeVelocity);
    plotAllTrials(fixEMobj, 3);
end

function plotAllTrials(fixEMobj, rowNo)
    
    nTrials = size(fixEMobj.emPosArcMin,1);
    visualizedSingleTrials = 5;
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', visualizedSingleTrials+1, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.01, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.01, ...
           'topMargin',      0.01);
       
    % Plot all trials
    xyRange = [-1 1]*max(abs(fixEMobj.emPosArcMin(:)));
    if (rowNo == 1)
        hFig = figure(1); clf;
        set(hFig, 'Position', [10 10 1600 900]);
    else
        figure(1);
    end
    
    
    for iTrial = 1:visualizedSingleTrials+1
        subplot('Position', subplotPosVectors(rowNo,iTrial).v);
        if (iTrial <= visualizedSingleTrials)
            plot(squeeze(fixEMobj.emPosArcMin(iTrial,:,1)), squeeze(fixEMobj.emPosArcMin(iTrial,:,2)), 'k-');
            title(sprintf('trial: %d', iTrial));
        else
            hold on
            for k = 1:nTrials
                plot(squeeze(fixEMobj.emPosArcMin(k,:,1)), squeeze(fixEMobj.emPosArcMin(k,:,2)), 'k-');
            end
            title(sprintf('%d trials\nsaccade type:\n''%s''', nTrials,fixEMobj.microSaccadeType));
        end
        hold on
        plot(xyRange, xyRange*0, 'r-');
        plot(xyRange*0, xyRange, 'r-');
        tickLabel = [xyRange(1) 0 xyRange(2)];
        set(gca, 'XLim', xyRange, 'YLim', xyRange, 'XTick', tickLabel, 'YTick', tickLabel, ...
                'XTickLabel', sprintf('%2.1f\n', tickLabel), 'YTickLabel', {});
        if (iTrial>1)
            set(gca, 'XTickLabel', {});
        else
            xlabel('arc min');
        end
        grid on
        axis 'square';
    end

end
