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
    set(ax.Legend, 'FontSize', ff.legendFontSize);
    set(ax.Legend, 'Box', ff.legendBox);

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
    
    drawnow;
end