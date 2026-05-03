function applyFormat(ax,ff)
    % Axis
    set(ax, ...
        'fontSize', ff.axisFontSize, ...
        'fontAngle', ff.axisFontAngle, ...
        'XColor', ff.axisColor, ...
        'YColor', ff.axisColor, ...
        'Color',ff.backgroundColor, ...
        'lineWidth', ff.axisLineWidth);

    % Title 
    set(ax.Title, ...
        'fontSize', ff.titleFontSize,...
        'fontWeight', ff.titleFontWeight, ...
        'color', ff.titleColor);

    % Legend
    set(ax.Legend, 'BackgroundAlpha', ff.legendBackgroundAlpha);
    set(ax.Legend, 'Color', ff.legendBackgroundColor);
    set(ax.Legend, 'TextColor', ff.legendTextColor);
    set(ax.Legend, 'EdgeColor', ff.legendEdgeColor);
    set(ax.Legend, 'LineWidth', ff.legendLineWidth);
    set(ax.Legend, 'FontSize', ff.legendFontSize);
    set(ax.Legend, 'Box', ff.legendBox);

    % Axes labels
    ax.XLabel.Color = ff.labelColor;
    ax.YLabel.Color = ff.labelColor;

    % Colorbar
    if (isfield(ff, 'colorbar'))
        set(ax.Colorbar, 'Location', ff.colorbar.location);
        set(ax.Colorbar, 'Orientation', ff.colorbar.orientation);
        set(ax.Colorbar, 'AxisLocation', ff.colorbar.axisLocation);
        set(ax.Colorbar, 'FontSize', ff.colorbar.fontSize);
        set(ax.Colorbar,'FontAngle', ff.colorbar.fontAngle);
    end

    % Ticks, grids, box
    xtickangle(ax, ff.axisTickAngle);
    set(ax, 'TickDir', ff.tickDir);
    grid(ax, ff.grid);
    set(ax, 'box', ff.box);
    

    try
        hFig = ancestor(ax, 'figure');
        set(hFig, 'Color', ff.backgroundColor)
    catch
        fprintf('PublicationReadeyPlotLib.applyFormat(): Unable to retrieve the parent figure handle from passed axes\n');
    end


    drawnow;
end