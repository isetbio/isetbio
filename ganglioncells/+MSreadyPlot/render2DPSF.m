function render2DPSF(ax, psfSupportXdegs, psfSupportYdegs, thePSFData, psfRangeDegs, plotTitle, ff, varargin)
    
    p = inputParser;
    p.addParameter('noXLabel', false, @islogical);
    p.addParameter('noYLabel', false, @islogical);
    p.addParameter('psfAlpha', 0.7, @(x)(isscalar(x)&&(x>0.0)&&(x<=1.0)));
    p.addParameter('withConeApertureData', [], @(x)(isempty(x)||isstruct(x)));
    p.addParameter('tickSeparationArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('gridlessPSF', false, @islogical);
    p.addParameter('superimposeXprofile', false, @islogical);
    p.addParameter('xProfileRange', [], @isscalar);
    p.addParameter('colorMap', [], @(x)(isempty(x)||(size(x,2) == 3)));
    p.parse(varargin{:});
    
    noXLabel = p.Results.noXLabel;
    noYLabel = p.Results.noYLabel;
    psfAlpha = p.Results.psfAlpha;
    coneApertureData = p.Results.withConeApertureData;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    gridlessPSF = p.Results.gridlessPSF;
    colorMap = p.Results.colorMap;
    superimposeXprofile = p.Results.superimposeXprofile;
    xProfileRange = p.Results.xProfileRange;
    psfRangeArcMin = psfRangeDegs*60;

    % plot
    psfSupportXarcmin = psfSupportXdegs * 60;
    psfSupportYarcmin = psfSupportYdegs * 60;
    if (isempty(colorMap))
        colorMap = brewermap(1024, 'blues');
    end


    if (~isempty(coneApertureData))
        xOutline = cosd(0:15:360);
        yOutline = sind(0:15:360);
        for iCone = 1:numel(coneApertureData.RcDegs)
            xx = coneApertureData.positionDegs(iCone,1)+coneApertureData.RcDegs(iCone)*xOutline;
            yy = coneApertureData.positionDegs(iCone,2)+coneApertureData.RcDegs(iCone)*yOutline;
            
            patch(ax, xx*60,yy*60, [0.8 0.8 0.8]);
            hold(ax, 'on');
            plot(ax, xx*60,yy*60, 'k-', 'LineWidth', 1.0, 'Color', [0.4 0.4 0.4]);
        end
        
    end

    zData = thePSFData/max(thePSFData(:));
    zLevels = [0.1 0.3 0.5 0.7 0.9];
    contourLineColor = [1 1 1]*0.3;
    cMosaic.semiTransparentContourPlot(ax, psfSupportXarcmin, psfSupportYarcmin, ...
        zData, zLevels, colorMap, psfAlpha, contourLineColor, 'LineWidth', 1.5)
    
    if (superimposeXprofile)
        hold(ax, 'on');
        [~, idx] = max(thePSFData(:));
        [mRow, ~] = ind2sub(size(thePSFData), idx);
        xProfile = squeeze(thePSFData(mRow,:));
        if (isempty(xProfileRange))
            xProfile = xProfile/max(xProfile(:));
        else
            xProfile = xProfile/xProfileRange;
        end

        plot(ax, psfSupportXarcmin, -psfRangeArcMin/2 + psfRangeArcMin * xProfile, 'k-', 'LineWidth', 1.5);
        hold(ax, 'off');
    end

    axis(ax, 'image');
    axis(ax, 'xy');

    
    if (isempty(tickSeparationArcMin))
        psfTicksArcMin = MSreadyPlot.spatialMapTicksArcMin(psfRangeArcMin);
    else
        psfTicksArcMin = -(10*tickSeparationArcMin):tickSeparationArcMin:(10*tickSeparationArcMin);
    end


    if (tickSeparationArcMin >= 6)
        psfTickLabels = sprintf('%2.0f\n', psfTicksArcMin);
    else
        psfTickLabels = sprintf('%2.1f\n', psfTicksArcMin);
    end

   
    % ticks and grids
    grid(ax, 'off'); box(ax, 'on');
    set(ax, 'XLim', 0.5*psfRangeArcMin*1.05*[-1 1], 'YLim', 0.5*psfRangeArcMin*1.05*[-1 1]);
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

    if (~gridlessPSF)
        grid(ax, 'on');
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
    colormap(ax, colorMap);
end
