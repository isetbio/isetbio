function removeRGCsOnPatchPerimeter(obj)
% Remove RGCs on the perimeter

    obj.visualizeCurrentConnectivityState(5555);

    allRGCindices = 1:size(obj.RGCRFcentroidsFromInputs,1);

    centerPos = obj.RGCRFcentroidsFromInputs;
    ss = squeeze(sum(obj.coneConnectivityMatrix,1));
    zeroInputRGCindices = find(ss == 0);
    centerPos(zeroInputRGCindices,:) = obj.RGCRFpositionsMicrons(zeroInputRGCindices,:);

    [~, RGCindicesToBeRemoved] = ...
        RGCconnector.pointsInsideBoundaryDefinedBySelectedPoints(centerPos, allRGCindices);

    currentRGCindices = 1:size(obj.coneConnectivityMatrix,2);
    remainingRGCindices = setdiff(currentRGCindices, RGCindicesToBeRemoved);
    fprintf('Removing %f RGCs on the perimeter\n', numel(RGCindicesToBeRemoved));

    % Disconnect cones from RGCs to be removed
    for iRGC = 1:numel(RGCindicesToBeRemoved)
        obj.coneConnectivityMatrix(:, RGCindicesToBeRemoved(iRGC)) = 0;
    end
    
    % Update connectivity matrix
    obj.coneConnectivityMatrix = obj.coneConnectivityMatrix(:, remainingRGCindices);


    obj.RGCRFpositionsMicrons = obj.RGCRFpositionsMicrons(remainingRGCindices,:);
    obj.RGCRFcentroidsFromInputs = obj.RGCRFcentroidsFromInputs(remainingRGCindices,:);
    obj.RGCRFspacingsMicrons = obj.RGCRFspacingsMicrons(remainingRGCindices);

    obj.visualizeCurrentConnectivityState(6666);
end