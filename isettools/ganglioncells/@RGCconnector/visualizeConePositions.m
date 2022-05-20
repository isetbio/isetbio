function visualizeConePositions(obj, ax, coneOutline)
% Visualize the cones of the input cone mosaic using a custom shape
% cone outline

    % Plot the cones
    coneTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID cMosaic.SCONE_ID];
    coneColors = [obj.inputConeMosaic.lConeColor; obj.inputConeMosaic.mConeColor; obj.inputConeMosaic.sConeColor];
    
    for iConeType = 1:numel(coneTypes)
        idx = find(obj.inputConeMosaic.coneTypes == coneTypes(iConeType));
        allConeRFpositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons(idx,:);
        allConeRFspacingsMicrons = obj.inputConeMosaic.coneRFspacingsMicrons(idx);
        [f,v] = RGCconnector.facesAndVertices(allConeRFpositionsMicrons, allConeRFspacingsMicrons, coneOutline);
        theColor = squeeze(coneColors(iConeType,:));
        patch(ax,'Faces', f, 'Vertices', v, 'FaceColor', theColor, 'EdgeColor', theColor*0.5);
        hold(ax, 'on')
    end
end
