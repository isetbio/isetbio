function visualizeInputPositions(obj)
% Visualize mRGC positions together with cone positions
%
% Syntax:
%   visualizeInputPositions(obj)
%
% Description:
%   Visualize mRGC positions together with imported ecc-varying cone positions
%   and together with cone positions in the equivalent employed reg-hex mosaic
%
% Inputs:
%    obj  - An instantiated @mRGCmosaic object
%
% Outputs:
%    none
%
% Optional key/value pairs:
%    none

% History:
%    11/06/2020  NPC   Wrote it
  
    % Compute xyRanges
    xRange = [...
         min(obj.importedData.rgcRFpositionsDegs(:,1))-0.05
         max(obj.importedData.rgcRFpositionsDegs(:,1))+0.05];
    yRange = [...
         min(obj.importedData.rgcRFpositionsDegs(:,2))-0.05
         max(obj.importedData.rgcRFpositionsDegs(:,2))+0.05];
     
    hFig = figure(); clf;
    set(hFig, 'Color', [1 1 1]);
    
    % Plot the cone positions in the ecc-varying cone mosaic together with the mRGC positions
    ax = subplot(1,2,1);
    plotData(ax, obj.importedData.coneRFpositionsDegs, obj.importedData.rgcRFpositionsDegs, ...
        xRange, yRange, 'imported cone positions')

    % Plot the cone positions in the regular hex cone mosaic together with the mRGC positions
    ax = subplot(1,2,2);
    plotData(ax, obj.inputConeMosaicMetaData.conePositionsDegs, obj.importedData.rgcRFpositionsDegs, ...
        xRange, yRange, 'regular hex mosaic cone positions')
end


function plotData(ax, coneRFpositionsDegs, rgcRFpositionsDegs, xRange, yRange, plotTitle)
    xOutline = cosd(0:15:360);
    yOutline = sind(0:15:360);
    
    coneFaceColor = [1 0.5 0.5];
    coneEdgeColor = [1 0.0 0.0];

    faceAlpha = 0.8;
    edgeAlpha = 0.7;
    
    % Compute radii of RFs based on their spacings, which are computed from their positions
    coneRadii = 0.48 * RGCmodels.Watson.convert.positionsToSpacings(coneRFpositionsDegs);
    rgcRadii = 0.48 * RGCmodels.Watson.convert.positionsToSpacings(rgcRFpositionsDegs);
    
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
