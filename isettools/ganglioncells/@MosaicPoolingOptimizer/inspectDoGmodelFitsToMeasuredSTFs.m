function inspectDoGmodelFitsToMeasuredSTFs(computeReadyMosaicFilename, mRGCMosaicSTFresponsesFilename)
    % Load the computed mRGC  STF responses 
    load(mRGCMosaicSTFresponsesFilename, ...
            'theMRGCMosaicOptimalSTFs', ...
            'visualRcDegsEstimates', ...
            'theMRGCMosaicSTFresponses', ...
            'theMRGCresponseTemporalSupportSeconds', ...
            'orientationsTested', 'spatialFrequenciesTested', ...
            'spatialPhasesDegs', 'coneContrasts');
    
    
    hFig = figure(1); clf;
    ff = MSreadyPlot.figureFormat('1x2 wide');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    axRhoSigma = theAxes{1,1};
    axRcDegs = theAxes{1,2};

    theSTFdata = theMRGCMosaicOptimalSTFs{1};
    idxKc = find(strcmp(theSTFdata.DoGfitParams.names, 'Kc'));
    idxKsKcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'kS/kC'));
    idxRsRcRatio = find(strcmp(theSTFdata.DoGfitParams.names, 'RsToRc'));
    idxRcDegs = find(strcmp(theSTFdata.DoGfitParams.names, 'RcDegs'));

    KsKcRatios = zeros(1, numel(theMRGCMosaicOptimalSTFs));
    RsRcRatios = zeros(1, numel(theMRGCMosaicOptimalSTFs));
    RcDegs = zeros(1, numel(theMRGCMosaicOptimalSTFs));

    for rgcIndex = 1:numel(theMRGCMosaicOptimalSTFs)
         theSTFdata = theMRGCMosaicOptimalSTFs{rgcIndex};  
         KsKcRatios(rgcIndex) = theSTFdata.DoGfitParams.finalValues(idxKsKcRatio);
         RsRcRatios(rgcIndex) = theSTFdata.DoGfitParams.finalValues(idxRsRcRatio);
         RcDegs(rgcIndex) = theSTFdata.DoGfitParams.finalValues(idxRcDegs);
    end
    % Integrated sensitivity ratios
    sigmaRatios = KsKcRatios .* (RsRcRatios).^2;

    scatter(axRhoSigma, RsRcRatios, sigmaRatios, 100, ...
        'MarkerFaceAlpha', 0.01, 'MarkerEdgeAlpha', 0.0, ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);
    xlabel(axRhoSigma, sprintf('\\rho'));
    ylabel(axRhoSigma, sprintf('\\sigma'));
    axis(axRhoSigma, 'square');
    set(axRhoSigma, 'XLim', [1 20], 'YLim', [0 1], 'YTick', 0:0.1:1, 'XTick', 0:2:20);
    set(axRhoSigma, 'FontSize', ff.fontSize, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    grid(axRhoSigma, 'on');
    box(axRhoSigma, 'off');

    scatter(axRcDegs, RcDegs*60, sigmaRatios, 100, ...
        'MarkerFaceAlpha', 0.01, 'MarkerEdgeAlpha', 0.0, ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 0 0]);

    xlabel(axRcDegs, 'Rc (arc min)');
    ylabel(axRcDegs, sprintf('\\sigma'));
    axis(axRcDegs, 'square');
    set(axRcDegs, 'XLim', [0 5], 'YLim', [0 1], 'YTick', 0:0.1:1, 'XTick', 0:1:20);
    set(axRcDegs, 'FontSize', ff.fontSize, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    grid(axRcDegs, 'on');
    box(axRcDegs, 'off');

    drawnow

    hFig = figure(2); clf;
    ff = MSreadyPlot.figureFormat('1x2 wide');
    theAxes = MSreadyPlot.generateAxes(hFig,ff);
    set(hFig, 'Color', [1 1 1]);
    axSTF = theAxes{1,1};
    axRMSE = theAxes{1,2};

    RMSE = nan(1, numel(theMRGCMosaicOptimalSTFs));
    maxRMSE = 0.0;
    
    

    for rgcIndex = 1:numel(theMRGCMosaicOptimalSTFs)
         theSTFdata = theMRGCMosaicOptimalSTFs{rgcIndex};

         % Compute RMS
         residuals = theSTFdata.DoGfit.compositeSTF(:) - theSTFdata.measured(:);
         RMSE(rgcIndex) = sqrt(mean(residuals.^2));
         edges = 0.01:0.002:0.10;
         h = histogram(axRMSE, RMSE, edges);
         hold(axRMSE, 'on');
         plot(axRMSE, maxRMSE*[1 1], [0 max(h.Values)], 'r-', 'LineWidth', ff.lineWidth);
         hold(axRMSE, 'off');
         axis(axRMSE, 'square')
         set(axRMSE, 'YTickLabel', {}, 'XLim', [0.01 0.1], 'XTick', 0.01:0.01:0.10, 'XTickLabel', {'','.02','','.04','','.06','','.08','','.1'}, ...
             'FontSize', ff.fontSize, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
         ylabel(axRMSE, {});
         xtickangle(axRMSE, 0);
         grid(axRMSE, 'on');
         box(axRMSE, 'off');
         xlabel(axRMSE, 'RMSE');
         title(axRMSE, sprintf('%d / %dmRGCs', rgcIndex, numel(theMRGCMosaicOptimalSTFs)));

         if (RMSE(rgcIndex) > maxRMSE)
             maxRMSE = RMSE(rgcIndex);
             cla(axSTF);
             MSreadyPlot.renderSTF(axSTF, spatialFrequenciesTested, theSTFdata.measured, ...
                 theSTFdata.DoGfit.sfHiRes, theSTFdata.DoGfit.compositeSTFHiRes, ...
                 theSTFdata.DoGfit.centerSTFHiRes, theSTFdata.DoGfit.surroundSTFHiRes, ...
                 sprintf('RGC %d (RMSE: %2.3f)', rgcIndex, maxRMSE), ...
                 [], ff, ...
                 'noYLabel', true);
         end

         drawnow;
         

    end

end
