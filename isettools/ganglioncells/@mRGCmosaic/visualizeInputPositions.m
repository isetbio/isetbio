function visualizeInputPositions(coneRFpositionsDegs, rgcRFpositionsDegs, coneRFpositionsDegsInRegHexMosaic)
          
    % Compute xyRanges
    xRange = [...
         min(rgcRFpositionsDegs(:,1))-0.05
         max(rgcRFpositionsDegs(:,1))+0.05];
    yRange = [...
         min(rgcRFpositionsDegs(:,2))-0.05
         max(rgcRFpositionsDegs(:,2))+0.05];
     
    % Compute spacing from positions
    coneSpacingHexRegMosaic = RGCmodels.Watson.convert.positionsToSpacings(coneRFpositionsDegsInRegHexMosaic);
    coneSpacing = RGCmodels.Watson.convert.positionsToSpacings(coneRFpositionsDegs);
    rgcSpacing = RGCmodels.Watson.convert.positionsToSpacings(rgcRFpositionsDegs);
    
    hFig = figure(1); clf;
    set(hFig, 'Color', [1 1 1]);
    
    ax = subplot(1,2,1);
    plotData(ax, coneRFpositionsDegs, rgcRFpositionsDegs, coneSpacing, rgcSpacing, xRange, yRange, 'imported cone positions')

    ax = subplot(1,2,2);
    plotData(ax, coneRFpositionsDegsInRegHexMosaic, rgcRFpositionsDegs, coneSpacingHexRegMosaic, rgcSpacing, xRange, yRange, 'regular hex mosaic cone positions')
end


function plotData(ax, coneRFpositionsDegs, rgcRFpositionsDegs, coneSpacings, rgcSpacings, xRange, yRange, plotTitle)
    xOutline = cosd(0:15:360);
    yOutline = sind(0:15:360);
    
    coneFaceColor = [1 0.5 0.5];
    coneEdgeColor = [1 0.0 0.0];

    faceAlpha = 0.8;
    edgeAlpha = 0.7;
    
    coneRadii = 0.48 * coneSpacings;
    rgcRadii = 0.48 * rgcSpacings;
    
    hold(ax, 'on');
    for k = 1:size(coneRFpositionsDegs,1)
        patchContour(ax, coneRFpositionsDegs(k,1) + coneRadii (k)*xOutline, ...
            coneRFpositionsDegs(k,2) + coneRadii (k)*yOutline, ...
            coneFaceColor, coneEdgeColor, faceAlpha, edgeAlpha);
    end
    for k = 1:size(rgcRFpositionsDegs,1)
        plot(ax, rgcRFpositionsDegs(k,1) + rgcRadii(k)*xOutline, ...
              rgcRFpositionsDegs(k,2) + rgcRadii(k)*yOutline, 'k-', ...
              'lineWidth', 2.0);
    end
    axis(ax, 'equal');
    set(ax, 'XLim', xRange, 'YLim', yRange);
    title(ax, plotTitle);
end

function patchContour(theAxes, xOutline, yOutline, faceColor, edgeColor, faceAlpha, edgeAlpha)
    v = [xOutline(:) yOutline(:)];
    f = 1:numel(xOutline);
    patch(theAxes, 'Faces', f, 'Vertices', v, 'FaceColor', faceColor, ...
            'FaceAlpha', faceAlpha, 'EdgeColor', edgeColor, ... 
           'EdgeAlpha', edgeAlpha, 'LineWidth', 1.5);
end
