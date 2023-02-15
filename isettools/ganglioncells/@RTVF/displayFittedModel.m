function displayFittedModel(figNo, modelConstants, theVisualRFmap, theVisualSTFdata, ...
    targetParams)

    visualRFcenterRcDegs = targetParams.RcDegs;
    targetRsRcatio = targetParams.RsRcRatio;
    targetSCIntSensRatio = targetParams.SCintSensRatio;
    
    figure(figNo); clf;
    ax = subplot(2,3,1);
    imagesc(ax,modelConstants.spatialSupportDegs(:,1), ...
            modelConstants.spatialSupportDegs(:,2), ...
            theVisualRFmap);
    axis(ax, 'image');  axis(ax, 'xy');

    title(sprintf('center max: %2.3f', max(theVisualRFmap(:))))

    ax = subplot(2,3,2);
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.visualSTF, 'ks', 'LineWidth', 1.5);
    hold(ax, 'on');
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.fittedDoGModelToVisualSTF.compositeSTF, 'r-', 'LineWidth', 1.5);
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.fittedDoGModelToVisualSTF.centerSTF, 'r:', 'LineWidth', 1.0);
    plot(ax,theVisualSTFdata.spatialFrequencySupport, theVisualSTFdata.fittedDoGModelToVisualSTF.surroundSTF, 'r--', 'LineWidth', 1.0);
    set(ax, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100]);

    ax = subplot(2,3,4);
    plot(ax,visualRFcenterRcDegs*60, theVisualSTFdata.fittedRcDegs*60 , 'ro', 'LineWidth', 1.5);
    hold(ax, 'on')
    plot([0 20], [0 20], 'k-', 'LineWidth', 1.0)
    xlabel('target Rc (arcmin)');
    ylabel('fitted DoG model Rc (arcmin)')
    set(ax, 'XLim', [0 3], 'YLim', [0 3]);
    axis(ax, 'square');

    ax = subplot(2,3,5);
    plot(ax,targetRsRcatio, theVisualSTFdata.fittedDoGModelRsRcRatio, 'ro', 'LineWidth', 1.5);
    hold(ax, 'on')
    plot([1 20], [1 20], 'k-', 'LineWidth', 1.0)
    xlabel('target Rs/Rc');
    ylabel('fitted DoG model Rs/Rc')
    set(ax, 'XLim', [1 20], 'YLim', [1 20]);
    axis(ax, 'square');

    ax = subplot(2,3,6);
    plot(ax, targetSCIntSensRatio, theVisualSTFdata.fittedDoGModelSCIntSensRatio , 'ro', 'LineWidth', 1.5);
    set(ax, 'XLim', [0 1], 'YLim', [0 1]);
    hold(ax, 'on')
    plot([0 1], [0 1], 'k-', 'LineWidth', 1.0)
    xlabel('target S/C int sens');
    ylabel('fitted S/C int sens')
    axis(ax, 'square');

end