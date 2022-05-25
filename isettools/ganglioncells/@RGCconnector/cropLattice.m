function cropLattice(obj, RGCRFposMicrons)
    
    % Crop the MRC mosaic so that we have enough space for the surround
    % cones
    maxRsRcRatio = 7;
    maxSurroundDiameter = maxRsRcRatio*prctile(obj.RGCRFspacingsMicrons,95);
    maxSurroundRadius = 0.5*maxSurroundDiameter;

    allConePositions = obj.inputConeMosaic.coneRFpositionsMicrons;
    minConePosX = min(allConePositions(:,1));
    minConePosY = min(allConePositions(:,2));
    maxConePosX = max(allConePositions(:,1));
    maxConePosY = max(allConePositions(:,2));

    idx = find(...
        (RGCRFposMicrons(:,1) >= minConePosX+maxSurroundRadius) & ...
        (RGCRFposMicrons(:,1) <= maxConePosX-maxSurroundRadius) & ...
        (RGCRFposMicrons(:,2) >= minConePosY+maxSurroundRadius) & ...
        (RGCRFposMicrons(:,2) <= maxConePosY-maxSurroundRadius));

    if (numel(idx) < 7)
        fprintf(2, 'Consider increasing the size of input cone mosaic\n');
        % Select the center most RGC
        d = sort(sum((bsxfun(@minus, RGCRFposMicrons, mean(allConePositions))).^2,2));
        idx = d(1:7);
    end

    % Crop positions
    obj.RGCRFpositionsMicrons = RGCRFposMicrons(idx,:);
       
    % Update spacings
    obj.RGCRFspacingsMicrons = obj.RGCRFspacingsMicrons(idx);

    % Initialize centroids
    obj.RGCRFcentroidsFromInputs = obj.RGCRFpositionsMicrons;
end