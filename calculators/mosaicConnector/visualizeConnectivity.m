function visualizeConnectivity(figNo, figureName, conePositionsMicrons, RGCRFPositionsMicrons, connectionMatrix, meanConesToRGCratio)
    hFig = figure(figNo); clf;
    set(hFig, 'Name', figureName);
    theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.03, ...
            'bottomMargin', 0.05, ...
            'rightMargin', 0.03, ...
            'topMargin', 0.03);
        
    % Cones in blue
    scatter(theAxesGrid{1,1}, conePositionsMicrons(:,1), conePositionsMicrons(:,2), 'b'); 
    for coneIndex = 1:size(conePositionsMicrons,1)
        text(conePositionsMicrons(coneIndex,1)+0.5, conePositionsMicrons(coneIndex,2)+1, sprintf('%d', coneIndex), 'Color', 'b');
    end
    hold(theAxesGrid{1,1}, 'on');
    % RGCs in black
    scatter(theAxesGrid{1,1}, RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 300,  [0.4 0.4 0.4]);
    for rgcIndex = 1:size(RGCRFPositionsMicrons,1)
        text(RGCRFPositionsMicrons(rgcIndex,1)+0.5, RGCRFPositionsMicrons(rgcIndex,2)+1, sprintf('%d', rgcIndex), 'Color', 'k');
    end
    
    rgcsNum = size(RGCRFPositionsMicrons,1);
    for theRGCindex = 1:rgcsNum  
        rgcPos = RGCRFPositionsMicrons(theRGCindex,:);

        % Find cones connected to this RGC
        coneWeights = squeeze(connectionMatrix(:, theRGCindex, 1));
        coneIndicesConnectedToThisRGC = find(coneWeights > 0);
        coneWeights = coneWeights(coneIndicesConnectedToThisRGC);
        lineColors = coneWeights / max(coneWeights);

        for iCone = 1:numel(coneIndicesConnectedToThisRGC)
            theConeIndex = coneIndicesConnectedToThisRGC(iCone);
            conePos = conePositionsMicrons(theConeIndex,:);
            line([rgcPos(1) conePos(1)], [rgcPos(2) conePos(2)], ...
                'LineWidth', 1.5, 'Color', (1-lineColors(iCone))*[1 1 1]);
        end
    end
    title(sprintf('mean cone-to-RGC ratio: %2.2f', meanConesToRGCratio));    
end