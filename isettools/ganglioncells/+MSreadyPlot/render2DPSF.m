function render2DPSF(ax, psfSupportXdegs, psfSupportYdegs, thePSFData, psfRangeDegs, plotTitle, ff, varargin)
    
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('withConeApertureData', [] , @(x)(isempty(x)||isstruct(x)));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    coneApertureData = p.Results.withConeApertureData;

    % plot
    psfSupportXarcmin = psfSupportXdegs * 60;
    psfSupportYarcmin = psfSupportYdegs * 60;
    cmap = brewermap(1024, 'greys');

    if (~isempty(coneApertureData))
        xOutline = cosd(0:15:360);
        yOutline = sind(0:15:360);
        for iCone = 1:numel(coneApertureData.RcDegs)
            xx = coneApertureData.positionDegs(iCone,1)+coneApertureData.RcDegs(iCone)*xOutline;
            yy = coneApertureData.positionDegs(iCone,2)+coneApertureData.RcDegs(iCone)*yOutline;
            plot(ax, xx*60,yy*60, 'k-', 'LineWidth', 1.0, 'Color', [0.4 0.4 0.4]);
            if (iCone == 1)
                hold(ax, 'on')
            end
        end
        
    end

    zData = thePSFData/max(thePSFData(:));
    zLevels = [0.1 0.3 0.5 0.7 0.9];
    alpha = 0.7;
    contourLineColor = [1 0 0];
    cMosaic.semiTransparentContourPlot(ax, psfSupportXarcmin, psfSupportYarcmin, ...
        zData, zLevels, cmap, alpha, contourLineColor, 'LineWidth', 1.5)
    

    axis(ax, 'image');
    axis(ax, 'xy');

    psfRangeArcMin = psfRangeDegs*60;
    psfTicksArcMin = MSreadyPlot.spatialMapTicksArcMin(psfRangeArcMin);

    if (psfRangeArcMin <= 2)
        psfTickLabels = sprintf('%2.1f\n', psfTicksArcMin);
    else
        psfTickLabels = sprintf('%2.0f\n', psfTicksArcMin);
    end

   
    % ticks and grids
    grid(ax, 'off'); box(ax, 'on');
    set(ax, 'XLim', psfRangeArcMin*1.05*[-1 1], 'YLim', psfRangeArcMin*1.05*[-1 1]);
    set(ax, 'XTick', psfTicksArcMin , 'YTick', psfTicksArcMin, ...
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
