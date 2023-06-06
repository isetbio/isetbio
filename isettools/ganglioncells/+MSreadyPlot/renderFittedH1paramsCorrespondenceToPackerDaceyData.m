function renderFittedH1paramsCorrespondenceToPackerDaceyData(ax, ...
  retinalConePoolingParams, ff)

    cla(ax, 'reset');

    idx = find(strcmp(retinalConePoolingParams.names,  'VnVwRatio'));
    fittedModel.NWvolumeRatio = retinalConePoolingParams.finalValues(idx);
    
    idx = find(strcmp(retinalConePoolingParams.names,  'RnRwRatio'));
    fittedModel.RnarrowToRwideRatio = retinalConePoolingParams.finalValues(idx);

    plot(ax, MosaicPoolingOptimizer.PackerDacey2002_H1params.NWvolumeRatios, ...
             MosaicPoolingOptimizer.PackerDacey2002_H1params.RnarrowToRwideRatios, 'kh', ...
             'MarkerSize', ff.markerSize, 'MarkerFaceColor', [0.8 0.8 0.8], 'LineWidth', ff.lineWidth);
    hold(ax, 'on')

    scatter(ax, fittedModel.NWvolumeRatio, fittedModel.RnarrowToRwideRatio, (ff.markerSize+2)^2,'o', ...
             'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0], ...
             'LineWidth', ff.lineWidth*1.5);

    
    set(ax, 'XLim', [1.21*ff.axisOffsetFactor 1.21], 'YLim', [0.4*ff.axisOffsetFactor 0.4], 'XTick', 0:0.2:1.2, 'YTick', 0:0.1:1, 'TickDir', 'both');

    axis(ax, 'square');
    grid(ax, 'on');
    box(ax, 'off');
    xtickangle(ax, 0);

    legend(ax, {'fitted H1 model params', 'H1 cell model params (Packer & Dacey)'}, ...
        'FontSize', ff.legendFontSize, 'box', 'off', 'Color', [0.8 0.8 0.7]);

    ylabel(ax,'V_{narrow}/V_{wide}', 'FontAngle', ff.axisFontAngle);
    xlabel(ax, 'R_{narrow}/R_{wide}', 'FontAngle', ff.axisFontAngle);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    title(ax, 'correspondence to Packer & Dacey', 'fontSize', ff.titleFontSize, ...
     'Color', ff.titleColor, 'FontWeight', ff.titleFontWeight); 

end
