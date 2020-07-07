function [RGCpositions, RGCpositionsNormalized] = determineRGCPositionsFromCenterInputs(theConeMosaic, patchEccentricityMicrons, centerWeights)
     % Retrieve cone positions (microns), cone spacings, and cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone positions: add the patch center so as to align with ecc-varying full mRGC mosaic
    conePositionsMicrons = bsxfun(@plus, cmStruct.coneLocsMicrons, patchEccentricityMicrons);
    
    rgcsNum = size(centerWeights,2);
    RGCpositions = zeros(rgcsNum,2);
    for iRGC = 1:rgcsNum
        weights = full(squeeze(centerWeights(:, iRGC)));
        centerIndices = find(weights>0);
        cp = conePositionsMicrons(centerIndices,:);
        wNet = sum(weights(centerIndices));
        cp(:,1) = cp(:,1) .* weights(centerIndices);
        cp(:,2) = cp(:,2) .* weights(centerIndices);
        RGCpositions(iRGC,:) = sum(cp,1)/wNet;
    end
    
    % Normalize to [0..1]
    mins = min(RGCpositions, [], 1);
    maxs = max(RGCpositions, [], 1);
    RGCpositionsNormalized(:,1) = (RGCpositions(:,1)-mins(1))/(maxs(1)-mins(1));
    RGCpositionsNormalized(:,2) = (RGCpositions(:,2)-mins(2))/(maxs(2)-mins(2));
end