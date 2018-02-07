function t_fixationalEyeMovementsTypes
% Examine eye movement paths produced by different values of the
% 'microSaccadeType' parameter.
%

% History
%   02/06/18  npc  Wrote it.
%   02/07/18  npc  Comments.

    close all;
    
    % Generate eye movement data for 20, 1 second long trials, with a
    % sample time of 1 millisecond.
    emDurationSeconds = 1.0; sampleTimeSeconds = 1/1000; nTrials = 256;
    
    % Do not compute velocity of eye movements
    computeVelocity = false;
    useParfor = true;
    
    % Random seed to be used in all eye movement compute() calls
    randomSeed = 1;
    
    % Visualize the 5 first trials
    visualizedSingleTrials = 5;
    
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
    tickLabel = [xyRange(1):10:xyRange(2)];
     
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
            xPos = squeeze(fixEMobj.emPosArcMin(:,:,1));
            yPos = squeeze(fixEMobj.emPosArcMin(:,:,2));
            fixationPos = xyRange(1):0.5:xyRange(2);
            midRow = find(fixationPos == 0);
            fixationMap = histcounts2(xPos(:),yPos(:),fixationPos,fixationPos);
            fixationMap = fixationMap / max(fixationMap(:));
            imagesc(fixationPos,fixationPos,fixationMap);
            hold on;
            position = fixationPos(1:end-1);
            slice = squeeze(fixationMap(midRow,:));
            plot(position, xyRange(1)+slice*xyRange(2)*0.9, '-', 'Color', [1 1 1], 'LineWidth', 1.5);
            slice = squeeze(fixationMap(:,midRow));
            plot(xyRange(2)-slice*xyRange(2)*0.9, position, '-', 'Color', [1 1 1], 'LineWidth', 1.5);
        end
        hold on
        plot(xyRange, xyRange*0, 'g-');
        plot(xyRange*0, xyRange, 'g-');
        set(gca, 'XLim', xyRange, 'YLim', xyRange, 'XTick', tickLabel, 'YTick', tickLabel, ...
                'XTickLabel', sprintf('%2.1f\n', tickLabel), 'YTickLabel', {});
        if (iTrial<visualizedSingleTrials+2)
            set(gca, 'XTickLabel', {});
        else
            xlabel('arc min');
        end
        grid on; box on; axis 'square'; axis 'xy'
    end
    colormap(jet(1024))
end
