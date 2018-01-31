function t_fixationalEyeMovements

    emDurationSeconds = 2; sampleTimeSeconds = 1/1000; nTrials = 1*12;
    
    % Initialize object
    fixEMobj = fixationalEM();
    fixEMobj.setDefaultParams();
    fixEMobj.microSaccadeType = 'none';
    
    tic
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials);
    totalTime = toc;
    fprintf('Total time: %2.1f seconds, %2.3f seconds/trial\n', totalTime, totalTime/size(fixEMobj.emPosArcMin,1));
    plotAllTrials(fixEMobj, 1);
    
    fixEMobj.setDefaultParams();
    fixEMobj.microSaccadeType = 'stats based';
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials);
    plotAllTrials(fixEMobj, 2);
    
    fixEMobj.setDefaultParams();
    fixEMobj.microSaccadeType = 'heatmap/fixation based';
    fixEMobj.compute(emDurationSeconds, sampleTimeSeconds, nTrials);
    plotAllTrials(fixEMobj, 3);
    
    % Plot last trial
    lastTrialEMposArcMin = squeeze(fixEMobj.emPosArcMin(end,:,:));
    figure(100); clf; hold on
    plot(fixEMobj.timeAxis, lastTrialEMposArcMin(:,1), 'r-'); hold on;
    plot(fixEMobj.timeAxis, lastTrialEMposArcMin(:,2), 'b-'); hold on;
    set(gca, 'YLim', 60*[-0.2 0.2]);
    ylabel('position (arc min)');
    drawnow;
end

function plotAllTrials(fixEMobj, figNo)
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 6, ...
           'colsNum', 12, ...
           'heightMargin',   0.02, ...
           'widthMargin',    0.02, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.04);
       
    % Plot all trials
    xyRange = [-1 1]*max(abs(fixEMobj.emPosArcMin(:)));
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1400 770]);
    for iTrial = 1:size(fixEMobj.emPosArcMin,1)
        i = floor((iTrial-1)/12)+1;
        j = mod(iTrial-1,12)+1;
        subplot('Position', subplotPosVectors(i,j).v);
        hold on
        if (i == 6) && (j == 12)
            plot(squeeze(fixEMobj.emPosArcMin(:,:,1)), squeeze(fixEMobj.emPosArcMin(:,:,2)), 'k.-');
            plot(xyRange, xyRange*0, 'r-');
            plot(xyRange*0, xyRange, 'r-');
            tickLabel = [xyRange(1) 0 xyRange(2)];
            set(gca, 'XLim', xyRange, 'YLim', xyRange, 'XTick', tickLabel, 'YTick', tickLabel, ...
                'XTickLabel', sprintf('%2.1f\n', tickLabel), 'YTickLabel', sprintf('%2.1f\n', tickLabel));
            xlabel('arc min');
            title(sprintf('saccade type:\n''%s''',fixEMobj.microSaccadeType));
        else
            plot(squeeze(fixEMobj.emPosArcMin(iTrial,:,1)), squeeze(fixEMobj.emPosArcMin(iTrial,:,2)), 'k.-');
            plot(xyRange, xyRange*0, 'r-');
            plot(xyRange*0, xyRange, 'r-');
            set(gca, 'XLim', xyRange, 'YLim', xyRange, 'XTick', [], 'YTick', []);
        end
        axis 'square';
    end
   
end
