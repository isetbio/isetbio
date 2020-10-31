function visualizeRetinalContrastProfiles(theRetinalLMScontrastSequence, gaborSpatialFrequencyCPD, LMScontrast, figExportsDir)
        
    plotlabOBJ = setupPlotLab(0, 5, 6);
    hFig = figure(100); clf;
    
    midRow = round(size(theRetinalLMScontrastSequence,2)/2);
    colsNum = size(theRetinalLMScontrastSequence,3);
    margin = round(colsNum/4);
    colsToUse = margin:colsNum-margin;

    for frame = 1:1 %size(theRetinalLMScontrastSequence,1)
        
        line(1:numel(colsToUse),squeeze(theRetinalLMScontrastSequence(frame,midRow,colsToUse,1)), 'Color', [1 0 0], 'LineWidth', 2); hold on;
        line(1:numel(colsToUse),squeeze(theRetinalLMScontrastSequence(frame,midRow,colsToUse,2)), 'Color', [0 0.7 0], 'LineWidth', 2); 
        line(1:numel(colsToUse),squeeze(theRetinalLMScontrastSequence(frame,midRow,colsToUse,3)), 'Color', [ 0 0 1], 'LineWidth', 2); 

        set(gca, 'YLim', [-0.11 0.11], 'XLim', [1 numel(colsToUse)], 'XTick', [], 'YTick', 0.1*[-1:0.5:1]);
        axis 'square'
        title(sprintf('%4.2f c/deg', gaborSpatialFrequencyCPD));
    end
    plotlab.offsetAxes(gca);
    drawnow;
    exportFig = true;
    if (exportFig)
       plotlabOBJ.exportFig(hFig, 'pdf', sprintf('RetinalContrastProfiles_%2.2fcpd_LMS_%0.2f_%0.2f_%0.2f', gaborSpatialFrequencyCPD, LMScontrast(1), LMScontrast(2), LMScontrast(3)), figExportsDir);
    end
    setupPlotLab(-1);
end

function plotlabOBJ = setupPlotLab(mode, figWidthInches, figHeightInches)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'in', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 8, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', figWidthInches, ...
                'figureHeightInches', figHeightInches);
    else
        pause(2.0);
        plotlab.resetAllDefaults();
    end
end 