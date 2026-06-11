function render2DVisualRF(ax, ...
        spatialSupportDegsX, spatialSupportDegsY, ...
        RFmap,  ...
        xyTicks, XLims, YLims, reverseXDir, plotTitle, noYLabel,  ff)

    rfSensitivityLUT = brewermap(1024, '*RdBu');
    backgroundColor = rfSensitivityLUT(512,:);
    imagesc(ax, spatialSupportDegsX, spatialSupportDegsY, RFmap);

    hold(ax, 'on');
    contour(ax, spatialSupportDegsX, spatialSupportDegsY, RFmap, ...
        exp(-1)*[0.99 1], 'LineWidth', 1.0, 'Color', [0 0 0]);
    contour(ax, spatialSupportDegsX, spatialSupportDegsY, RFmap, ...
        -0.01*[0.99 1], 'LineWidth', 1.0, 'LineStyle', ':', 'Color', [0 0 0]);

    axis(ax, 'image');
    axis(ax, 'xy');
    set(ax, 'Color', backgroundColor, 'CLim', [-1 1], 'XLim', XLims, 'YLim', YLims, ...
        'XTick', xyTicks, 'YTick', xyTicks, ...
        'XTickLabel', sprintf('%2.1f\n', xyTicks), ...
        'YTickLabel', sprintf('%2.1f\n', xyTicks));
    colormap(ax,rfSensitivityLUT);
    xtickangle(ax, 0);
    
    if (reverseXDir)
        set(ax, 'XDir', 'reverse');
    end

    % Finalize subplot
    if (isempty(ff))
        xlabel(ax, 'space, x (deg)');
        if (~noYLabel)
            ylabel(ax, 'space, y (deg)');
        end
    else
        xlabel(ax, 'space, x (deg)', 'FontAngle', ff.axisFontAngle);
        if (~noYLabel)
            ylabel(ax, 'space, y (deg)', 'FontAngle', ff.axisFontAngle);
        end

        if (~isempty(plotTitle))
        title(ax, 'surround boosted (x10)', ...
            'FontSize', ff.titleFontSize, ...
            'FontWeight', ff.titleFontWeight, ...
            'Color', ff.titleColor);
        end

        % Font size
        set(ax, 'FontSize', ff.fontSize);

        % axis color and width
        set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end
end