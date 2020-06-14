function visualizeConnectivity(figNo, figureName, conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, coneTypes, meanConesToRGCratio)
    hFig = figure(figNo); clf;
    set(hFig, 'Name', figureName);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.01, ...
            'topMargin', 0.05);
        
    % Display cones
    LconeIndices = find(coneTypes == 2);
    MconeIndices = find(coneTypes == 3);
    SconeIndices = find(coneTypes == 4);
    hold(theAxesGrid{1,1}, 'on');
    scatter(theAxesGrid{1,1}, conePositionsMicrons(LconeIndices,1), conePositionsMicrons(LconeIndices,2), 'MarkerFaceColor', [1 0.5 0.5], 'MarkerEdgeColor', [1 0 0]);
    scatter(theAxesGrid{1,1}, conePositionsMicrons(MconeIndices,1), conePositionsMicrons(MconeIndices,2), 'MarkerFaceColor', [0.3 0.8 0.3], 'MarkerEdgeColor', [0 1 0]);
    scatter(theAxesGrid{1,1}, conePositionsMicrons(SconeIndices,1), conePositionsMicrons(SconeIndices,2), 'MarkerFaceColor', [0.5 0.5 1.0], 'MarkerEdgeColor', [0 0 1]);
    
    
    for coneIndex = 1:size(conePositionsMicrons,1)
        text(theAxesGrid{1,1}, conePositionsMicrons(coneIndex,1)+0.5, conePositionsMicrons(coneIndex,2)+1, sprintf('%d', coneIndex), 'Color', 'b');
    end
    
    % RGCs in black
    scatter(theAxesGrid{1,1}, RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 300,  [0.4 0.4 0.4]);
    for rgcIndex = 1:size(RGCRFPositionsMicrons,1)
        text(theAxesGrid{1,1}, RGCRFPositionsMicrons(rgcIndex,1)+0.5, RGCRFPositionsMicrons(rgcIndex,2)+1, sprintf('%d', rgcIndex), 'Color', 'k');
    end
    
    rgcsNum = size(RGCRFPositionsMicrons,1);
    for theRGCindex = 1:rgcsNum  
        rgcPos = RGCRFPositionsMicrons(theRGCindex,:);

        % Find cones connected to this RGC
        coneWeights = squeeze(connectionMatrix(:, theRGCindex));
        coneIndicesConnectedToThisRGC = find(coneWeights > 0);
        coneWeights = coneWeights(coneIndicesConnectedToThisRGC);
        lineColors = coneWeights / max(coneWeights);

        for iCone = 1:numel(coneIndicesConnectedToThisRGC)
            theConeIndex = coneIndicesConnectedToThisRGC(iCone);
            conePos = conePositionsMicrons(theConeIndex,:);
            line(theAxesGrid{1,1},[rgcPos(1) conePos(1)], [rgcPos(2) conePos(2)], ...
                'LineWidth', 1.5, 'Color', (1-lineColors(iCone))*[1 1 1]);
        end
    end
    title(theAxesGrid{1,1}, sprintf('mean cone-to-RGC ratio: %2.2f', meanConesToRGCratio));    
end