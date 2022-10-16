function generateCenterSurroundSpatialPoolingRF(obj, theRetinaToVisualFieldTransformerOBJList, ...
            theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
            theVisualSTFSurroundToCenterRcRatioGrid, ...
            theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid)

    obj.theRetinaToVisualFieldTransformerOBJList = theRetinaToVisualFieldTransformerOBJList;
    obj.theOpticsPositionGrid = theOpticsPositionGrid;
    obj.theConesNumPooledByTheRFcenterGrid = theConesNumPooledByTheRFcenterGrid;
    obj.theVisualSTFSurroundToCenterRcRatioGrid = theVisualSTFSurroundToCenterRcRatioGrid;
    obj.theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid = theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid;


    % Reset the center and surround cone pooling matrices (sparse)
    conesNum = size(obj.rgcRFcenterConeConnectivityMatrix,1);
    rgcsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    obj.rgcRFcenterConePoolingMatrix = sparse(conesNum, rgcsNum);
    obj.rgcRFsurroundConePoolingMatrix = sparse(conesNum, rgcsNum);
  
    for iRGC = 1:rgcsNum
        

        % Retrieve this cell's center cone indices and weights
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        indicesOfCenterCones = find(abs(connectivityVector) > 0.0001);

        % Compute distances of this RGC to all eccGrids
        distancesToEccGrid = sum((bsxfun(@minus, obj.theOpticsPositionGrid, obj.rgcRFpositionsDegs(iRGC,:))).^2,2);

        % Retrieve the correct RTVFTobj based on this cells position and
        % #of center cones. For now only checking the centerConesNum
%         iObj = find(...
%             (obj.theConesNumPooledByTheRFcenterGrid == numel(indicesOfCenterCones)) & ...  % match the conesNum in the center
%             (distancesToEccGrid == min(distancesToEccGrid)) ...                     % match to the closest eccGrid position
%             )

        iObj = find(...
            (obj.theConesNumPooledByTheRFcenterGrid == numel(indicesOfCenterCones)) ...  % match the conesNum in the center
            );
        theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};

        fprintf('Computing c/s cone pooling weights for RGC %d of %d using the %d RTVFTobj\n', iRGC, rgcsNum, iObj);
        
        % Extract the retinal cone pooling weights compute function
        theRetinalConePoolingWeightsComputeFunction = theRTVFTobj.rfComputeStruct.modelConstants.weightsComputeFunctionHandle;
    
        % Extract the retinal cone pooling weights params vector
        theRetinalConePoolingParamsVector = theRTVFTobj.rfComputeStruct.retinalConePoolingParams.finalValues;

        % Extract the model constants
        modelConstants = theRTVFTobj.rfComputeStruct.modelConstants;
        
        % Update the center cone indices and weights for the target midgetRGC center 
        modelConstants.indicesOfCenterCones = indicesOfCenterCones;
        modelConstants.weightsOfCenterCones = connectivityVector(indicesOfCenterCones);
    
        % Compute the pooled cone indices and weights using the current
        % theRetinalConePoolingWeightsComputeFunction &
        % theRetinalConePoolingParamsVector
        pooledConeIndicesAndWeights = theRetinalConePoolingWeightsComputeFunction(...
            modelConstants, theRetinalConePoolingParamsVector);
    
        % Update theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
        % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices
        centerIndices = pooledConeIndicesAndWeights.centerConeIndices;
        centerWeights = pooledConeIndicesAndWeights.centerConeWeights;
        surroundIndices = pooledConeIndicesAndWeights.surroundConeIndices;
        surroundWeights = pooledConeIndicesAndWeights.surroundConeWeights;
        obj.rgcRFcenterConePoolingMatrix(centerIndices,iRGC) = centerWeights;
        obj.rgcRFsurroundConePoolingMatrix(surroundIndices,iRGC) = surroundWeights;
    end


end