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
        obj.theRetinaToVisualFieldTransformerOBJList{1}.LconeRFcomputeStruct.modelConstants.coneWeightsCompensateForVariationsInConeEfficiency;

    if (compensateForVariationsInConeEfficiency)
        totalCenterStrengthForLconeCenterModelRF = zeros(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
        totalSurroundStrengthForLconeCenterModelRF = zeros(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
        totalCenterStrengthForMconeCenterModelRF = zeros(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
        totalSurroundStrengthForMconeCenterModelRF = zeros(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));

        for iObj = 1:numel(obj.theRetinaToVisualFieldTransformerOBJList)
            theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};

            % The L-cone center models
            idx = find(theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights > 0);
            weights = theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights(idx);
            totalCenterStrengthForLconeCenterModelRF(iObj) = sum(weights(:));
            idx = find(theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights > 0);
            weights = theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights(idx);
            totalSurroundStrengthForLconeCenterModelRF(iObj) = sum(weights(:));

            % The M-cone center models
            idx = find(theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights > 0);
            weights = theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights(idx);
            totalCenterStrengthForMconeCenterModelRF(iObj) = sum(weights(:));
            idx = find(theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights > 0);
            weights = theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights(idx);
            totalSurroundStrengthForMconeCenterModelRF(iObj) = sum(weights(:));
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

        % Determine whether the majority of center cones are L- or M-
        theMajorityCenterConeType = obj.majorityCenterConeType(iRGC);
        
        % Compute the indices of the triangulating RTVFobjects and their
        % contributing weights
        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
            obj.triangulatingRTVFobjectIndicesAndWeights(iRGC);

        % Compute the pooledConeIndicesAndWeights for each of the triangulatingRTVFobjIndices
        pooledConeIndicesAndWeightsForNearbyObj = cell(1, numel(triangulatingRTVFobjIndices));
        for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)

            iObj = triangulatingRTVFobjIndices(iNearbyObj);
            theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};

            switch (theMajorityCenterConeType)
                case cMosaic.LCONE_ID
                    totalCenterStrengthForModelRF = totalCenterStrengthForLconeCenterModelRF(iObj);
                    totalSurroundStrengthForModelRF = totalSurroundStrengthForLconeCenterModelRF(iObj);
                    theRFcomputeStruct = theRTVFTobj.LconeRFcomputeStruct;

                case cMosaic.MCONE_ID
                    totalCenterStrengthForModelRF = totalCenterStrengthForMconeCenterModelRF(iObj);
                    totalSurroundStrengthForModelRF = totalSurroundStrengthForMconeCenterModelRF(iObj);
                    theRFcomputeStruct = theRTVFTobj.MconeRFcomputeStruct;

                otherwise
                    error('How can the majority cone type be not L- or M- ??')
            end

            pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj} = computePooledConeIndicesAndWeightsFromRTVFTobject(...
                theRFcomputeStruct, indicesOfCenterCones, weightsOfCenterCones, ...
                totalCenterStrengthForModelRF, totalSurroundStrengthForModelRF, ...
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

function pooledConeIndicesAndWeights = computePooledConeIndicesAndWeightsFromRTVFTobject(theRFcomputeStruct, ...
    indicesOfCenterCones, weightsOfCenterCones,  ...
    totalCenterStrengthForModelRF, totalSurroundStrengthForModelRF, compensateForVariationsInConeEfficiency)

    % Extract the retinal cone pooling weights compute function
    theRetinalConePoolingWeightsComputeFunction = theRFcomputeStruct.modelConstants.weightsComputeFunctionHandle;
    
    % Extract the retinal cone pooling weights params vector
    theRetinalConePoolingParamsVector = theRFcomputeStruct.retinalConePoolingParams.finalValues;

    % Extract the model constants
    modelConstants = theRFcomputeStruct.modelConstants;
       
    % Extract the model center and surround retinal cone pooling weights
    fittedCenterConeWeights = theRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights;
    fittedSurroundConeWeights = theRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights;
    fittedNonConnectableSurroundConeWeights = theRFcomputeStruct.pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights;

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