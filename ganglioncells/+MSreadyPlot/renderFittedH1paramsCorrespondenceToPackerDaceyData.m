function renderFittedH1paramsCorrespondenceToPackerDaceyData(ax, ...
  retinalConePoolingParams, ff, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('resetAxes', true, @islogical);
    p.addParameter('noLegends', false, @islogical);
    p.parse(varargin{:});
    
    noLegends = p.Results.noLegends;

    if (p.Results.resetAxes)
        cla(ax, 'reset');
    end

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

    if (~noLegends)
        legend(ax, {'Packer & Dacey', 'fitted H1 model'}, ...
            'FontSize', ff.legendFontSize, 'box', 'off', 'Color', [0.8 0.8 0.7]);
    end

    ylabel(ax, 'Vnarrow / Vwide', 'FontAngle', ff.axisFontAngle);
    xlabel(ax, 'Rnarrow / Rwide', 'FontAngle', ff.axisFontAngle);

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);

end
