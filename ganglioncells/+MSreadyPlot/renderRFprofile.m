function renderRFprofile(ax, spatialSupport, theRFprofile, ...
        xyTicksArcMin, xyLimsArcMin, reverseXDir, gridlessLineWeightingFuncions, ...
        axisLabel, ff)

    shadedAreaPlot(ax,spatialSupport, theRFprofile, 0, [1 0.85 0.85], [0.7 0 0], 0.5, 1.5);
    hold(ax, 'on');
    plot(ax, spatialSupport, theRFprofile*0, 'k-', 'LineWidth', 1.0);

    set(ax, 'YLim', [-0.4 1], 'YTick', -1:0.2:1, ...
        'XTick', xyTicksArcMin, 'XLim', xyLimsArcMin, ...
        'XTickLabel', sprintf('%2.0f\n', xyTicksArcMin), ...
        'YTickLabel', {});
    axis(ax, 'square');

    if (gridlessLineWeightingFuncions)
        grid(ax, 'off');
    else
        grid(ax, 'on');
    end
    box(ax, 'on');
    
    xtickangle(ax, 0);

    if (reverseXDir)
        set(ax, 'XDir', 'reverse');
    end

    % Finalize subplot
    if (isempty(ff))
        xlabel(ax, axisLabel);
    else
        xlabel(ax, axisLabel, 'FontAngle', ff.axisFontAngle);

        % Font size
        set(ax, 'FontSize', ff.fontSize);

        % axis color and width
        set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    end
end



function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth)
    x = [x fliplr(x)];
    y = [y y*0+baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, 'FaceAlpha', faceAlpha, 'LineWidth', lineWidth);
end
