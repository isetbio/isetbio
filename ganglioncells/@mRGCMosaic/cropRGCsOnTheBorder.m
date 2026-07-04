function cropRGCsOnTheBorder(obj)
    % Remove RGC RFs on the margins 
    indicesOfRGCRFsToKeep = findRGCRFsWithinBorders(obj.rgcRFpositionsMicrons, obj.rgcRFspacingsMicrons);
    
    % Update rf center connectivity matrix & rgcRFpositionsMicrons
    obj.rgcRFcenterConeConnectivityMatrix = obj.rgcRFcenterConeConnectivityMatrix(:,indicesOfRGCRFsToKeep);

    obj.rgcRFpositionsMicrons = obj.rgcRFpositionsMicrons(indicesOfRGCRFsToKeep,:);
    obj.rgcRFspacingsMicrons = obj.rgcRFspacingsMicrons(indicesOfRGCRFsToKeep);
    obj.rgcRFpositionsDegs = obj.rgcRFpositionsDegs(indicesOfRGCRFsToKeep,:);
    obj.rgcRFspacingsDegs = obj.rgcRFspacingsDegs(indicesOfRGCRFsToKeep);

end

function indicesOfRFsToKeep = findRGCRFsWithinBorders(rgcRFpositions, rgcRFspacings)
    [minXY, idx] = min(rgcRFpositions,[],1);
    maxSpacing = max(rgcRFspacings(idx));
    [maxXY, idx] = max(rgcRFpositions,[],1);
    maxSpacing = max([maxSpacing max(rgcRFspacings(idx))]);

    xLims(1) = minXY(1)+maxSpacing;
    xLims(2) = maxXY(1)-maxSpacing;
    yLims(1) = minXY(2)+maxSpacing;
    yLims(2) = maxXY(2)-maxSpacing;
    indicesOfRFsToKeep= find(...
            (rgcRFpositions(:,1)>=xLims(1)) & ...
            (rgcRFpositions(:,1)<=xLims(2)) & ...
            (rgcRFpositions(:,2)>=yLims(1)) & ...
            (rgcRFpositions(:,2)<=yLims(2)) ...
            );
    
end
