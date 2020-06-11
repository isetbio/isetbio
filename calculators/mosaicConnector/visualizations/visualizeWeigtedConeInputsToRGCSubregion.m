function visualizeWeigtedConeInputsToRGCSubregion(theAxes, coneTypes, ...
    conePositionsMicrons, coneSpacingsMicrons, connectionWeights,  maxWeightVisualized, ...
    rgcPositionMicrons, rgcSubregionRadiusMicrons,  spatialSupportRangeMicrons)

    xx = cosd(0:5:360);
    yy = sind(0:5:360);
    
    connectionWeights = connectionWeights / maxWeightVisualized;
    connectionWeights = min([ones(numel(connectionWeights),1) connectionWeights],[],2);
    
    backgroundColor = [0.9 0.9 0.9];
    
    for iCone = 1:numel(connectionWeights)
        switch (coneTypes(iCone))
            case 2
               coneColor = [1 0 0]*connectionWeights(iCone) + backgroundColor * (1-connectionWeights(iCone));
            case 3
               coneColor = [0 1 0]*connectionWeights(iCone) + backgroundColor * (1-connectionWeights(iCone));
            case 4
               coneColor = [0 0 1]*connectionWeights(iCone) + backgroundColor * (1-connectionWeights(iCone));
        end
        coneColor = min([ones(1,3); coneColor]);
        cPosMicrons = conePositionsMicrons(iCone,:);
        cRadiusMicrons = 0.8*0.5*coneSpacingsMicrons(iCone);
        cProfile(:,1) = cPosMicrons(1) + cRadiusMicrons*xx;
        cProfile(:,2) = cPosMicrons(2) + cRadiusMicrons*yy;
        patch(theAxes, 'faces', 1:size(cProfile,1), 'Vertices', cProfile, ...
                'FaceColor', coneColor, 'FaceAlpha', 1, 'EdgeColor', [0.4 0.4 0.4], 'EdgeAlpha', 0.4);
        hold(theAxes, 'on');
    end % iCone
    
    % Plot the outline of the subregion
    rgcSubregionOutline(:,1) = rgcPositionMicrons(1) + rgcSubregionRadiusMicrons * xx;
    rgcSubregionOutline(:,2) = rgcPositionMicrons(2) + rgcSubregionRadiusMicrons * yy;

    patch(theAxes,'faces', 1:size(rgcSubregionOutline,1), 'Vertices', rgcSubregionOutline, ...
            'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.0, 'EdgeColor', [0 0 0], 'EdgeAlpha', 0.4);
        
    set(theAxes, 'XLim', rgcPositionMicrons(1) + spatialSupportRangeMicrons*[-1 1], ...
                 'YLim', rgcPositionMicrons(2) + spatialSupportRangeMicrons*[-1 1], ...
                 'Color', backgroundColor);
    
    axis(theAxes,'square');     
end