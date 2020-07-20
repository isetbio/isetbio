function  visualizeResponseComponentsForTargetRGC(targetRGC, responseTimeAxis, centerResponses, surroundResponses, ...
    spatialFrequenciesCPD, maxSpikeRate, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, figExportsDir)

    plotlabOBJ = setupPlotLab(0, 16, 12);
    hFig = figure(123); clf;
    
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 4, ...
        'colsNum', 4, ...
        'leftMargin', 0.06, ...
        'widthMargin', 0.03, ...
        'heightMargin', 0.03, ...
        'bottomMargin', 0.05, ...
        'rightMargin', 0.0, ...
        'topMargin', 0.01);
    
    centerResponsesMean = squeeze(mean(centerResponses,2));
    surroundResponsesMean = squeeze(mean(surroundResponses,2));
    
    for sfIndex = 1:numel(spatialFrequenciesCPD)  
        row = floor((sfIndex-1)/4)+1;
        col = mod(sfIndex-1,4)+1;
        ax = theAxesGrid{row,col};
        centerResponses = squeeze(centerResponsesMean(sfIndex,targetRGC,:));
        surroundResponses = squeeze(surroundResponsesMean(sfIndex,targetRGC,:));
        line(ax, responseTimeAxis, centerResponses, 'Color', [1 0 0], 'LineWidth', 1.5); hold on;
        line(ax, responseTimeAxis, surroundResponses, 'Color', [0 0 1], 'LineWidth', 1.5);
        set(ax, 'XLim', [0 0.5], 'YLim', maxSpikeRate*[0 1], 'XTick', 0:0.1:0.5, 'YTick',(0:0.2:1)*maxSpikeRate);
        if (row ==3)&&(col == 1)
            xlabel(ax, 'time (sec)');
            ylabel(ax, 'response');
        else
            set(ax, 'XTickLabel', {}, 'YTickLabel', {});
        end
        axis(ax, 'square');
        text(ax,0.02, 180,sprintf('%4.2f c/deg', spatialFrequenciesCPD(sfIndex)), 'FontSize',16);
        
    end
    
    drawnow;
    if (exportFig)
       pdfFileName = sprintf('ResponseComponents_RGC_%2.0f_LMS_%0.2f_%0.2f_%0.2f_PolansSID_%d_%s', targetRGC, LMScontrast(1), LMScontrast(2), LMScontrast(3), PolansSubjectID, opticsPostFix);
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