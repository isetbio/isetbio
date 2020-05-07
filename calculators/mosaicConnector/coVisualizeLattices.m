function coVisualizeLattices(conePositionsMicrons, RGCRFPositionsMicrons)
    
    maxEcc = 5;
    minX = 2; maxX = 4;
    eccRangeY = 0.6;
    
    micronsPerDeg = 300;
    ecc = sqrt(sum(conePositionsMicrons.^2,2));
    idx = find(ecc<1.1*maxEcc*micronsPerDeg);
    conePositionsMicrons = conePositionsMicrons(idx,:);
    
    ecc = sqrt(sum(RGCRFPositionsMicrons.^2,2));
    idx = find(ecc<1.1*maxEcc*micronsPerDeg);
    RGCRFPositionsMicrons = RGCRFPositionsMicrons(idx,:);
    
    % Instantiate a plotlab object
    plotlabOBJ = plotlab();
    
    figWidth = 30;
    figHeight = figWidth * eccRangeY / (maxX-minX);
    plotlabOBJ.applyRecipe(...
            'colorOrder', [1 0 0; 0 0 0], ...
            'figureWidthInches', figWidth, ...
            'figureHeightInches', figHeight);
        
    hFig = figure(1); clf;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
            'leftMargin', 0.05, ...
            'bottomMargin', 0.05, ...
            'rightMargin', 0.07, ...
            'topMargin', 0.05);
    theAxesGrid = theAxesGrid{1,1};
    hold(theAxesGrid, 'on');
    
    visualizeLattice(conePositionsMicrons);
    visualizeLattice(RGCRFPositionsMicrons);
    axis 'equal';
    set(gca, 'XLim', [minX maxX]*micronsPerDeg, ...
        'YLim', 0.5*eccRangeY*[-1 1]*micronsPerDeg);
    
    
end

function visualizeLattice(rfPositions)
    x = rfPositions(:,1);
    y = rfPositions(:,2);
    scatter(x,y,100);
    
    showLattice = ~true;
    if (showLattice)
        triangleIndices = delaunayn(rfPositions);

        xx = []; yy = [];
        for triangleIndex = 1:size(triangleIndices, 1)
            coneIndices = triangleIndices(triangleIndex, :);
            xCoords = x(coneIndices);
            yCoords = y(coneIndices);
            for k = 1:numel(coneIndices)
                xx = cat(2, xx, xCoords);
                yy = cat(2, yy, yCoords);
            end
        end

        patch(xx, yy, [0 0 1], 'EdgeColor', color, ...
            'EdgeAlpha', 0.8, 'FaceAlpha', 0.0, ...
            'FaceColor', [0.8 0.8 0.8], 'LineWidth', 1.0, ...
            'LineStyle', '-', 'Parent', gca); 
        hold on;
    end
    
end
