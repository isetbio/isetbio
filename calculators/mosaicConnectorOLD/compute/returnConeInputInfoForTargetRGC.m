function targetRGCconeInputInfo = returnConeInputInfoForTargetRGC(targetRGCindices, theMidgetRGCmosaic, theConeMosaic)
    
    % Retrieve all cone types
    cmStruct = theConeMosaic.geometryStructAlignedWithSerializedConeMosaicResponse();
    
    % Cone types
    coneTypes = cmStruct.coneTypes;
    
    targetRGCconeInputInfo = cell(1, numel(targetRGCindices));
    
    for mRGCindex = 1:numel(targetRGCindices)
        targetRGCindex = targetRGCindices(mRGCindex);
        % Retrieve indices of cones connected to the RF center
        centerWeights = full(squeeze(theMidgetRGCmosaic.centerWeights(:, targetRGCindex)));
        centerConeIndices = find(centerWeights > 0);

        % Retrieve indices of cones connected to the RF surround
        surroundWeights = full(squeeze(theMidgetRGCmosaic.surroundWeights(:, targetRGCindex)));
        surroundConeIndices = find(surroundWeights > 0);

        % Return types, weights and indices of cones connected to the RF center
        coneInputInfo.mRGCindex = targetRGCindex;
        coneInputInfo.center.types = coneTypes(centerConeIndices);
        coneInputInfo.center.weights = centerWeights(centerConeIndices);
        coneInputInfo.center.indices = centerConeIndices;

        % Return types, weights and indices of cones connected to the RF surround
        coneInputInfo.surround.types = coneTypes(surroundConeIndices);
        coneInputInfo.surround.weights = surroundWeights(surroundConeIndices);
        coneInputInfo.surround.indices = surroundConeIndices;
        
        targetRGCconeInputInfo{mRGCindex} = coneInputInfo;
    end
    
end

