function generateCenterSurroundSpatialPoolingRFs(obj, theRetinaToVisualFieldTransformerOBJList, ...
            theSamplingPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
            theVisualSTFSurroundToCenterRcRatioGrid, ...
            theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid)

    obj.theRetinaToVisualFieldTransformerOBJList = theRetinaToVisualFieldTransformerOBJList;
    obj.theSamplingPositionGrid = theSamplingPositionGrid;
    obj.theConesNumPooledByTheRFcenterGrid = theConesNumPooledByTheRFcenterGrid;
    obj.theVisualSTFSurroundToCenterRcRatioGrid = theVisualSTFSurroundToCenterRcRatioGrid;
    obj.theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid = theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid;


    % Reset the center and surround cone pooling (sparse) matrices
    conesNum = size(obj.rgcRFcenterConeConnectivityMatrix,1);
    rgcsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    obj.rgcRFcenterConePoolingMatrix = sparse(conesNum, rgcsNum);
    obj.rgcRFsurroundConePoolingMatrix = sparse(conesNum, rgcsNum);

    % Was the model center/surround weights fitted with compensation for
    % varions in efficacy with eccentricity?
    compensateForVariationsInConeEfficiency = ...
        obj.theRetinaToVisualFieldTransformerOBJList{1}.rfComputeStruct.modelConstants.coneWeightsCompensateForVariationsInConeEfficiency;

    if (compensateForVariationsInConeEfficiency)
        totalCenterStrengthForModelRF = zeros(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
        totalSurroundStrengthForModelRF = zeros(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));

        for iObj = 1:numel(obj.theRetinaToVisualFieldTransformerOBJList)
            theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};
            idx = find(theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeWeights > 0);
            weights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeWeights(idx);
            totalCenterStrengthForModelRF(iObj) = sum(weights(:));
           
            idx = find(theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeWeights > 0);
            weights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeWeights(idx);
            totalSurroundStrengthForModelRF(iObj) = sum(weights(:));
        end
    end

    % Sort RGCs with respect to their distance from the mRGC mosaic center
    distancesToEccGrid = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mean(obj.rgcRFpositionsDegs,1))).^2,2);
    [~,sortedRGCindices] = sort(distancesToEccGrid, 'ascend');

    for iSortedRGCindex = 1:numel(sortedRGCindices)
        fprintf('Generating surround for RGC%d of %d RGC\n', iSortedRGCindex, numel(sortedRGCindices));

        % Retrieve RGC index
        iRGC = sortedRGCindices(iSortedRGCindex);

        % Retrieve this cell's # of center cone indices
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        indicesOfCenterCones = find(connectivityVector > 0.0001);
        weightsOfCenterCones = connectivityVector(indicesOfCenterCones);

        % Compute the indices of the triangulating RTVFobjects and their
        % contributing weights
        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
            obj.triangulatingRTVFobjectIndicesAndWeights(iRGC);

        % Compute the pooledConeIndicesAndWeights for each of the triangulatingRTVFobjIndices
        pooledConeIndicesAndWeightsForNearbyObj = cell(1, numel(triangulatingRTVFobjIndices));
        for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)
            iObj = triangulatingRTVFobjIndices(iNearbyObj);
            theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};
            pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj} = computePooledConeIndicesAndWeightsFromRTVFTobject(...
                theRTVFTobj, indicesOfCenterCones, weightsOfCenterCones, ...
                totalCenterStrengthForModelRF(iObj), totalSurroundStrengthForModelRF(iObj), ...
                compensateForVariationsInConeEfficiency);
        end
        
        % Accumulate weigted surround weights
        for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)
            newSurroundConeIndices = pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.surroundConeIndices;
            newSurroundConeWeights = pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.surroundConeWeights * triangulatingRTVFobjWeights(iNearbyObj);

            if (iNearbyObj == 1) 
                surroundConeIndices = newSurroundConeIndices;
                surroundConeWeights = newSurroundConeWeights;
            else
                % Find which of the newSurroundConeIndices already exist in
                % the surroundConeIndices
                [ia,ib] = ismember(newSurroundConeIndices, surroundConeIndices);

                % Sum weights (to the weights of previously existing cone indices)
                for i = 1:numel(ia)
                    if (ia(i) == 1)
                        surroundConeWeights(ib(i)) = surroundConeWeights(ib(i)) + newSurroundConeWeights(i); 
                    end
                end
                % Add weights (of not previoulsy included cone indices)
                for i = 1:numel(ia)
                    if (ia(i) == 0)
                        surroundConeWeights = cat(1, surroundConeWeights, newSurroundConeWeights(i));
                        surroundConeIndices = cat(1, surroundConeIndices, newSurroundConeIndices(i));
                    end
                end
            end
        end
        
        % Pooled cone indices and weights for this RGC
        pooledConeIndicesAndWeights.centerConeWeights = pooledConeIndicesAndWeightsForNearbyObj{1}.centerConeWeights;
        pooledConeIndicesAndWeights.centerConeIndices = pooledConeIndicesAndWeightsForNearbyObj{1}.centerConeIndices;

        pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
        pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;

        % Set the theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
        % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices
        obj.rgcRFcenterConePoolingMatrix(pooledConeIndicesAndWeights.centerConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.centerConeWeights;
        obj.rgcRFsurroundConePoolingMatrix(pooledConeIndicesAndWeights.surroundConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.surroundConeWeights;
      
        if (max(pooledConeIndicesAndWeights.centerConeWeights(:)) ~= 1)
            error('Max center weight must be 1');
        end

    end
end

function pooledConeIndicesAndWeights = computePooledConeIndicesAndWeightsFromRTVFTobject(theRTVFTobj, ...
    indicesOfCenterCones, weightsOfCenterCones,  ...
    totalCenterStrengthForModelRF, totalSurroundStrengthForModelRF, compensateForVariationsInConeEfficiency)

    % Extract the retinal cone pooling weights compute function
    theRetinalConePoolingWeightsComputeFunction = theRTVFTobj.rfComputeStruct.modelConstants.weightsComputeFunctionHandle;
    
    % Extract the retinal cone pooling weights params vector
    theRetinalConePoolingParamsVector = theRTVFTobj.rfComputeStruct.retinalConePoolingParams.finalValues;

    % Extract the model constants
    modelConstants = theRTVFTobj.rfComputeStruct.modelConstants;
       
    % Extract the model center and surround retinal cone pooling weights
    fittedCenterConeWeights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeWeights;
    fittedSurroundConeWeights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeWeights;
    fittedNonConnectableSurroundConeWeights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights;

    % Update the center cone indices and weights for the target midgetRGC center 
    % so we can pass that info to theRetinalConePoolingWeightsComputeFunction()
    modelConstants.indicesOfCenterCones = indicesOfCenterCones;
    modelConstants.weightsOfCenterCones = weightsOfCenterCones;


    % Compute the pooled cone indices and weights using the
    % theRetinalConePoolingWeightsComputeFunction & theRetinalConePoolingParamsVector
    pooledConeIndicesAndWeights = theRetinalConePoolingWeightsComputeFunction(...
            modelConstants, theRetinalConePoolingParamsVector);

    % Extract the model center and surround retinal cone pooling weights
    currentRGCCenterConeWeights = pooledConeIndicesAndWeights.centerConeWeights;
    currentRGCSurroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights;
    currentRGCNonConnectableSurroundConeWeights = pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights;


    if (compensateForVariationsInConeEfficiency)
        % Adjust for total strength of the center
        totalCenterStrengthForThisRF = sum(pooledConeIndicesAndWeights.centerConeWeights);
        centerStrengthCorrectionFactor = totalCenterStrengthForModelRF/totalCenterStrengthForThisRF;
        pooledConeIndicesAndWeights.centerConeWeights = pooledConeIndicesAndWeights.centerConeWeights * centerStrengthCorrectionFactor;

        % Surround strength correction factor
        % Adjust for different # of S-cones
        fittedBoost = 1 + sum(fittedNonConnectableSurroundConeWeights)/sum(fittedSurroundConeWeights);
        currentBoost = 1 + sum(currentRGCNonConnectableSurroundConeWeights)/sum(currentRGCSurroundConeWeights);
        surroundBoostRatio = currentBoost/fittedBoost;
        pooledConeIndicesAndWeights.surroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights * surroundBoostRatio;

        % Adjust for total strength of the surround
        totalSurroundStrengthForThisRF = sum(pooledConeIndicesAndWeights.surroundConeWeights);
        surroundStrengthCorrectionFactor = totalSurroundStrengthForModelRF/totalSurroundStrengthForThisRF;
        pooledConeIndicesAndWeights.surroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights * surroundStrengthCorrectionFactor;
     end % compensateForVariationsInConeEfficiency

     % Normalize weights
     normalizingWeight = max(pooledConeIndicesAndWeights.centerConeWeights(:));
     pooledConeIndicesAndWeights.centerConeWeights = pooledConeIndicesAndWeights.centerConeWeights / normalizingWeight;
     pooledConeIndicesAndWeights.surroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights / normalizingWeight;

end