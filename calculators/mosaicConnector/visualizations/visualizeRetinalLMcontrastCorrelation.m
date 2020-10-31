function visualizeRetinalLMcontrastCorrelation(spatialFrequenciesCPD, roiLcontrast, roiMcontrast, ...
    LMScontrast, opticsPostFix, PolansSubjectID, exportFig, figExportsDir)
    plotlabOBJ = setupPlotLab(0, 16, 12);
    hFig = figure(101); clf;
    
    rowsNum = 3;
    colsNum = 5;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', rowsNum, ...
        'colsNum', colsNum, ...
        'leftMargin', 0.06, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.03, ...
        'bottomMargin', 0.05, ...
        'rightMargin', 0.0, ...
        'topMargin', 0.01);
    
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        row = floor((sfIndex-1)/colsNum)+1;
        col = mod(sfIndex-1,colsNum)+1;
        ax = theAxesGrid{row,col};
        retinalLconeContrast = roiLcontrast{sfIndex};
        retinalMconeContrast = roiMcontrast{sfIndex};
        scatter(ax, retinalLconeContrast(1:4:end), retinalMconeContrast(1:4:end), '.');
        set(ax, 'XLim', 0.1*[-1 1], 'YLim', 0.1*[-1 1], 'XTick', 0.1*(-1:0.5:1), 'YTick', 0.1*[-1:0.5:1]);
        if (row ==3)&&(col == 1)
            xlabel(ax, 'retinal L-cone contrast');
            ylabel(ax, 'retinal M-cone contrast');
        else
            set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        end
        axis(ax, 'square');
        text(ax,-0.09, 0.085,sprintf('%4.2f c/deg', spatialFrequenciesCPD(sfIndex)), 'FontSize',16);
        drawnow;
    end
    drawnow;
    if (exportFig)
        pdfFileName = sprintf('RetinalLMContrastCorrelationAtDifferentSFs_LMS_%0.2f_%0.2f_%0.2f_PolansSID_%d_%s', LMScontrast(1), LMScontrast(2), LMScontrast(3), PolansSubjectID, opticsPostFix);
        
        plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, figExportsDir);
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