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

        % Determine the net L-cone, and net M-cone weights in the RF center of this RGC
        theCenterConeTypeWeights = obj.centerConeTypeWeights(iRGC);
        theTotalCenterConeWeights = sum(theCenterConeTypeWeights([cMosaic.LCONE_ID cMosaic.MCONE_ID]));

        % Weights for L- and M-cone compute structs proportional to the
        % # of center L- and M-cones, respectively
        theLconeComputeStructWeight = theCenterConeTypeWeights(cMosaic.LCONE_ID) / theTotalCenterConeWeights;
        theMconeComputeStructWeight = theCenterConeTypeWeights(cMosaic.MCONE_ID) / theTotalCenterConeWeights;

        % Compute the indices of the triangulating RTVFobjects and their
        % contributing weights
        [triangulatingRTVFobjIndices, triangulatingRTVFobjWeights] = ...
            obj.triangulatingRTVFobjectIndicesAndWeights(iRGC);


        % Compute the pooledConeIndicesAndWeights for each of the triangulatingRTVFobjIndices
        pooledLConeIndicesAndWeightsForNearbyObj = cell(1, numel(triangulatingRTVFobjIndices));
        pooledMConeIndicesAndWeightsForNearbyObj = cell(1, numel(triangulatingRTVFobjIndices));

        for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)

            iObj = triangulatingRTVFobjIndices(iNearbyObj);
            theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};

            if (theLconeComputeStructWeight > 0) 
                pooledLConeIndicesAndWeightsForNearbyObj{iNearbyObj} = computePooledConeIndicesAndWeightsFromRTVFTobject(...
                    theRTVFTobj.LconeRFcomputeStruct, indicesOfCenterCones, weightsOfCenterCones, ...
                    totalCenterStrengthForLconeCenterModelRF(iObj), totalSurroundStrengthForLconeCenterModelRF(iObj), ...
                    compensateForVariationsInConeEfficiency);
            end

            if (theMconeComputeStructWeight > 0) 
                pooledMConeIndicesAndWeightsForNearbyObj{iNearbyObj} = computePooledConeIndicesAndWeightsFromRTVFTobject(...
                    theRTVFTobj.MconeRFcomputeStruct, indicesOfCenterCones, weightsOfCenterCones, ...
                    totalCenterStrengthForMconeCenterModelRF(iObj), totalSurroundStrengthForMconeCenterModelRF(iObj), ...
                    compensateForVariationsInConeEfficiency);
            end
        end


        % Accumulate weigted center weights from Lcone and Mcone compute
        % struct contributions
        centerConeWeights = zeros(size(indicesOfCenterCones));

        for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)
            if (theLconeComputeStructWeight > 0) 
                centerConeWeights = centerConeWeights + ...
                    pooledLConeIndicesAndWeightsForNearbyObj{iNearbyObj}.centerConeWeights * triangulatingRTVFobjWeights(iNearbyObj) * theLconeComputeStructWeight;
            end
            if (theMconeComputeStructWeight > 0) 
                centerConeWeights = centerConeWeights + ...
                    pooledMConeIndicesAndWeightsForNearbyObj{iNearbyObj}.centerConeWeights * triangulatingRTVFobjWeights(iNearbyObj) * theMconeComputeStructWeight;
            end
        end


        % Initialize surround cone indices and weights
        surroundConeIndices = [];
        surroundConeWeights = [];

        % Accumulate weighted surround weights: Lcone compute struct contribution
        [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
            theLconeComputeStructWeight, ...
            pooledLConeIndicesAndWeightsForNearbyObj, ...
            triangulatingRTVFobjWeights, ...
            surroundConeIndices, surroundConeWeights);


        % Accumulate weighted surround weights: Mcone compute struct contribution
        [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
            theMconeComputeStructWeight, ...
            pooledMConeIndicesAndWeightsForNearbyObj, ...
            triangulatingRTVFobjWeights, ...
            surroundConeIndices, surroundConeWeights);


        % Pooled cone indices and weights for this RGC        
        pooledConeIndicesAndWeights.centerConeIndices = indicesOfCenterCones;
        pooledConeIndicesAndWeights.centerConeWeights = centerConeWeights;

        pooledConeIndicesAndWeights.surroundConeIndices = surroundConeIndices;
        pooledConeIndicesAndWeights.surroundConeWeights = surroundConeWeights;

        % Set the theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
        % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices
        obj.rgcRFcenterConePoolingMatrix(pooledConeIndicesAndWeights.centerConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.centerConeWeights;
        obj.rgcRFsurroundConePoolingMatrix(pooledConeIndicesAndWeights.surroundConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.surroundConeWeights;
      
        if (abs((max(pooledConeIndicesAndWeights.centerConeWeights(:)) - 1)) > 0.001)
            error('Max center weight must be 1. It is %2.4f', max(pooledConeIndicesAndWeights.centerConeWeights(:)));
        end

    end
end

function [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
            theComputeStructWeight, ...
            pooledConeIndicesAndWeightsForNearbyObj, ...
            triangulatingRTVFobjWeights, ...
            surroundConeIndices, surroundConeWeights)

    % Accumulate weigted surround weights: Lcone compute struct contribution
    if (theComputeStructWeight > 0) 
        for iNearbyObj = 1:numel(pooledConeIndicesAndWeightsForNearbyObj)
            newSurroundConeIndices = pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.surroundConeIndices;
            newSurroundConeWeights = pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.surroundConeWeights * triangulatingRTVFobjWeights(iNearbyObj) * theComputeStructWeight;
    
            % Find which of the newSurroundConeIndices already exist in the surroundConeIndices
            [ia,ib] = ismember(newSurroundConeIndices, surroundConeIndices);
    
            for i = 1:numel(ia)
                if (ia(i) == 0)
                    % Include weights for previoulsy NOT included cone indices
                    surroundConeWeights = cat(1, surroundConeWeights, newSurroundConeWeights(i));
                    surroundConeIndices = cat(1, surroundConeIndices, newSurroundConeIndices(i));
                else
                    % Accumulate weights for previously included cone indices
                    surroundConeWeights(ib(i)) = surroundConeWeights(ib(i)) + newSurroundConeWeights(i); 
                end
            end
                
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