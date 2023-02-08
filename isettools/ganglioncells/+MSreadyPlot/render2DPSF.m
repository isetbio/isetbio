function render2DPSF(ax, psfSupportXdegs, psfSupportYdegs, thePSFData, psfRangeDegs, plotTitle, ff, varargin)
    
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;

    % plot
    psfSupportXarcmin = psfSupportXdegs * 60;
    psfSupportYarcmin = psfSupportYdegs * 60;
    cmap = brewermap(1024, 'greys');
    MSreadyPlot.transparentContourPlot(ax, ...
        psfSupportXarcmin , psfSupportYarcmin , ...
        thePSFData/max(thePSFData(:)), ...
        [0.1 0.3 0.5 0.7 0.9], cmap);
    axis(ax, 'image');

    psfTicksMin = (-30:5:30);
    psfRangeArcMin = psfRangeDegs*60;
    if (psfRangeArcMin <= 2)
        psfTicks = (-3:0.5:3);
    elseif (psfRangeArcMin <= 5)
        psfTicks = 0.1*psfTicksMin;
    elseif (psfRangeArcMin <= 10)
        psfTicks = 0.4*psfTicksMin;
    elseif (psfRangeArcMin <= 20)
        psfTicks = 1*psfTicksMin;
    elseif (psfRangeArcMin <= 40)
        psfTicks = 2*psfTicksMin;
    elseif (psfRangeArcMin <= 50)
        psfTicks = 3*psfTicksMin; 
    elseif (psfRangeArcMin <= 60)
        psfTicks = 4*psfTicksMin; 
    elseif (psfRangeArcMin <= 100)
        psfTicks = 6*psfTicksMin; 
    elseif (psfRangeArcMin <= 200)
        psfTicks = 10*psfTicksMin;
    elseif (psfRangeArcMin <= 400)
        psfTicks = 20*psfTicksMin; 
    end

    if (psfRangeArcMin <= 2)
        psfTickLabels = sprintf('%2.1f\n', psfTicks);
    else
        psfTickLabels = sprintf('%2.0f\n', psfTicks);
    end

   
    % ticks and grids
    grid(ax, 'on'); box(ax, 'off');
    set(ax, 'XLim', psfRangeArcMin*1.05*[-1 1], 'YLim', psfRangeArcMin*1.05*[-1 1]);
    set(ax, 'XTick', psfTicks , 'YTick', psfTicks, ...
        'XTickLabel', psfTickLabels, 'YTickLabel', psfTickLabels);
    set(ax, 'TickDir', 'in')

    % xy axis labels
    xtickangle(ax, 0);
    ytickangle(ax, 0);
    if (~noXLabel)
        xlabel(ax, 'space, x (arc min)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'XTickLabel', {});
    end
    if (~noYLabel)
        ylabel(ax, 'space, y (arc min)', 'FontAngle', ff.axisFontAngle);
    else
        set(ax, 'YTickLabel', {});
    end

    % Font size
    set(ax, 'FontSize', ff.fontSize);

    % axis color and width
    set(ax, 'XColor', ff.axisColor, 'YColor', ff.axisColor, 'LineWidth', ff.axisLineWidth);
    
    % plot title
    if (~isempty(plotTitle))
        title(ax, plotTitle, 'Color', ff.titleColor, 'FontSize', ff.titleFontSize, 'FontWeight', ff.titleFontWeight);
    end

    % Colormap
    colormap(ax, cmap);

end
