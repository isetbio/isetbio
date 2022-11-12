function wireRFcentersToInputCones(obj, visualizeAlignment, visualizeWiringStages)

    % Extract positions and spacings of cones lying within the FOV of the RGC mosaic
    % Note that there are cones located outside the FOV of the RGC mosaic which feed to
    % the surrounds of the peripheral RGCs
    [coneRFPositionsDegs, coneRFPositionsMicrons, ...
     coneRFSpacingsDegs, coneRFSpacingsMicrons, idxConesInsideRGCmosaic] = conesWithinTheRGCmosaic(obj);
 
    
    % Step1. Align each RGC with its nearest cone. 
    % Alignment occurs only for RGCs for which the cone-to-RGC ratio < 2. 
    % This ensures that all RGC's in central retina are connected to at least one cone. 
    % Since cones are more numerous than RGCs, some cones will not connect to an RGC at this step. 
    obj.alignRGCsWithCones(coneRFPositionsMicrons, coneRFPositionsDegs, visualizeAlignment);

    % Step 2. Connect L and M cones to RGC centers
    obj.connectConesToRGCcenters(idxConesInsideRGCmosaic, visualizeWiringStages);
    
end

function [coneRFpositionsDegs, coneRFpositionsMicrons, ...
     coneRFspacingsDegs, coneRFspacingsMicrons, idxConesToBeConnected] = conesWithinTheRGCmosaic(obj)

    % Indices of cones to be connected to the centers of RGCs
    idxConesToBeConnected = find(...
        (obj.inputConeMosaic.coneRFpositionsMicrons(:,1) >= obj.minRFpositionMicrons(1)) & ...
        (obj.inputConeMosaic.coneRFpositionsMicrons(:,2) >= obj.minRFpositionMicrons(2)) & ...
        (obj.inputConeMosaic.coneRFpositionsMicrons(:,1) <= obj.maxRFpositionMicrons(1)) & ...
        (obj.inputConeMosaic.coneRFpositionsMicrons(:,2) <= obj.maxRFpositionMicrons(2)));
    
    % Indices of cones that are not to be connected to the RGC centers
    conesNum = size(obj.inputConeMosaic.coneRFpositionsMicrons,1);
    idxConesOutsideRGCmosaic = setdiff(1:conesNum, idxConesToBeConnected);
    
    % Ensure that these cones do not get connected to the RGC centers.
    % We do this by setting their position/spacing to Inf
    n = numel(idxConesOutsideRGCmosaic);
    coneRFpositionsMicrons = obj.inputConeMosaic.coneRFpositionsMicrons;
    coneRFpositionsDegs = obj.inputConeMosaic.coneRFpositionsDegs;
    coneRFspacingsDegs = obj.inputConeMosaic.coneRFspacingsDegs;
    coneRFspacingsMicrons = obj.inputConeMosaic.coneRFspacingsMicrons;
    
    coneRFpositionsDegs(idxConesOutsideRGCmosaic,:) = Inf(n,2);
    coneRFpositionsMicrons(idxConesOutsideRGCmosaic,:) = Inf(n,2);
    coneRFspacingsMicrons(idxConesOutsideRGCmosaic) = Inf(n,1);
    coneRFspacingsDegs(idxConesOutsideRGCmosaic) = Inf(n,1);
end

