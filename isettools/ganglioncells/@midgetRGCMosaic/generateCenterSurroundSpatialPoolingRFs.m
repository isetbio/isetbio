function generateCenterSurroundSpatialPoolingRFs(obj, theRTVFobjList , ...
            theSamplingPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
            theVisualSTFSurroundToCenterRcRatioGrid, ...
            theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid)

    obj.theRTVFobjList = theRTVFobjList ;
    obj.theSamplingPositionGrid = theSamplingPositionGrid;
    obj.theConesNumPooledByTheRFcenterGrid = theConesNumPooledByTheRFcenterGrid;
    obj.theVisualSTFSurroundToCenterRcRatioGrid = theVisualSTFSurroundToCenterRcRatioGrid;
    obj.theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid = theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid;


    % Reset the center and surround cone pooling (sparse) matrices
    conesNum = size(obj.rgcRFcenterConeConnectivityMatrix,1);
    rgcsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    obj.rgcRFcenterConePoolingMatrix = sparse(conesNum, rgcsNum);
    obj.rgcRFsurroundConePoolingMatrix = sparse(conesNum, rgcsNum);

    compensateForVariationsInConeEfficiency = false;
    if (compensateForVariationsInConeEfficiency)
        totalCenterStrengthForLconeCenterModelRF = zeros(1, numel(obj.theRTVFobjList ));
        totalSurroundStrengthForLconeCenterModelRF = zeros(1, numel(obj.theRTVFobjList ));
        totalCenterStrengthForMconeCenterModelRF = zeros(1, numel(obj.theRTVFobjList ));
        totalSurroundStrengthForMconeCenterModelRF = zeros(1, numel(obj.theRTVFobjList ));
    
        for iObj = 1:numel(obj.theRTVFobjList )
            theRTVFTobj = obj.theRTVFobjList {iObj};
    
            % Center and surround strengths for the L-cone center model
            idx = find(theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights > 0);
            weights = theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights(idx);
            totalCenterStrengthForLconeCenterModelRF(iObj) = sum(weights(:));
            idx = find(theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights > 0);
            weights = theRTVFTobj.LconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights(idx);
            totalSurroundStrengthForLconeCenterModelRF(iObj) = sum(weights(:));
    
             % Center and surround strengths for the M-cone center model
            idx = find(theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights > 0);
            weights = theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.centerConeWeights(idx);
            totalCenterStrengthForMconeCenterModelRF(iObj) = sum(weights(:));
            idx = find(theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights > 0);
            weights = theRTVFTobj.MconeRFcomputeStruct.pooledConeIndicesAndWeights.surroundConeWeights(idx);
            totalSurroundStrengthForMconeCenterModelRF(iObj) = sum(weights(:));
        end
    else
        totalCenterStrengthForLconeCenterModelRF = nan(1, numel(obj.theRTVFobjList ));
        totalSurroundStrengthForLconeCenterModelRF = nan(1, numel(obj.theRTVFobjList ));
        totalCenterStrengthForMconeCenterModelRF = nan(1, numel(obj.theRTVFobjList ));
        totalSurroundStrengthForMconeCenterModelRF = nan(1, numel(obj.theRTVFobjList ));
    end

    

    % Sort RGCs with respect to their distance from the mRGC mosaic center
    distancesToEccGrid = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mean(obj.rgcRFpositionsDegs,1))).^2,2);
    [~,sortedRGCindices] = sort(distancesToEccGrid, 'ascend');

    centerConeIndicesAllRGCs = cell(rgcsNum,1);
    centerConeWeightsAllRGCs = cell(rgcsNum,1);
    centerRGCindicesAllRGCs = cell(rgcsNum,1);

    surroundConeIndicesAllRGCs = cell(rgcsNum,1);
    surroundConeWeightsAllRGCs = cell(rgcsNum,1);
    surroundRGCindicesAllRGCs = cell(rgcsNum,1);

    parfor iSortedRGCindex = 1:numel(sortedRGCindices)
        fprintf('Generating surround for RGC%d of %d RGC\n', iSortedRGCindex, numel(sortedRGCindices));

        % Retrieve RGC index
        iRGC = sortedRGCindices(iSortedRGCindex);

        % Retrieve this cell's # of center cone indices
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        centerConeIndicesForThisRGC = find(connectivityVector > 0.0001);
        weightsOfCenterConesForThisRGC = connectivityVector(centerConeIndicesForThisRGC);

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
            theRTVFTobj = obj.theRTVFobjList{iObj};

            if (theLconeComputeStructWeight > 0) 
                pooledLConeIndicesAndWeightsForNearbyObj{iNearbyObj} = computePooledConeIndicesAndWeightsFromRTVFTobject(...
                    theRTVFTobj.LconeRFcomputeStruct, centerConeIndicesForThisRGC , weightsOfCenterConesForThisRGC, ...
                    totalCenterStrengthForLconeCenterModelRF(iObj), totalSurroundStrengthForLconeCenterModelRF(iObj));
            end

            if (theMconeComputeStructWeight > 0) 
                pooledMConeIndicesAndWeightsForNearbyObj{iNearbyObj} = computePooledConeIndicesAndWeightsFromRTVFTobject(...
                    theRTVFTobj.MconeRFcomputeStruct, centerConeIndicesForThisRGC, weightsOfCenterConesForThisRGC, ...
                    totalCenterStrengthForMconeCenterModelRF(iObj), totalSurroundStrengthForMconeCenterModelRF(iObj));
            end
        end


        % Accumulate weigted center weights from Lcone and Mcone compute
        % struct contributions
        centerConeWeightsForThisRGC = zeros(size(centerConeIndicesForThisRGC));

        for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)
            if (theLconeComputeStructWeight > 0) 
                centerConeWeightsForThisRGC = centerConeWeightsForThisRGC + ...
                    pooledLConeIndicesAndWeightsForNearbyObj{iNearbyObj}.centerConeWeights * triangulatingRTVFobjWeights(iNearbyObj) * theLconeComputeStructWeight;
            end
            if (theMconeComputeStructWeight > 0) 
                centerConeWeightsForThisRGC = centerConeWeightsForThisRGC + ...
                    pooledMConeIndicesAndWeightsForNearbyObj{iNearbyObj}.centerConeWeights * triangulatingRTVFobjWeights(iNearbyObj) * theMconeComputeStructWeight;
            end
        end


        % Initialize surround cone indices and weights
        surroundConeIndicesForThisRGC = [];
        surroundConeWeightsForThisRGC = [];

        % Accumulate weighted surround weights: Lcone compute struct contribution
        [surroundConeIndicesForThisRGC, surroundConeWeightsForThisRGC] = accumulateSurroundWeights(...
            theLconeComputeStructWeight, ...
            pooledLConeIndicesAndWeightsForNearbyObj, ...
            triangulatingRTVFobjWeights, ...
            surroundConeIndicesForThisRGC, surroundConeWeightsForThisRGC);


        % Accumulate weighted surround weights: Mcone compute struct contribution
        [surroundConeIndicesForThisRGC, surroundConeWeightsForThisRGC] = accumulateSurroundWeights(...
            theMconeComputeStructWeight, ...
            pooledMConeIndicesAndWeightsForNearbyObj, ...
            triangulatingRTVFobjWeights, ...
            surroundConeIndicesForThisRGC, surroundConeWeightsForThisRGC);


        % Pooled cone indices and weights for this RGC 
        pooledConeIndicesAndWeightsForThisRGC = struct();
        pooledConeIndicesAndWeightsForThisRGC.centerConeIndices = centerConeIndicesForThisRGC;
        pooledConeIndicesAndWeightsForThisRGC.centerConeWeights = centerConeWeightsForThisRGC;

        pooledConeIndicesAndWeightsForThisRGC.surroundConeIndices = surroundConeIndicesForThisRGC;
        pooledConeIndicesAndWeightsForThisRGC.surroundConeWeights = surroundConeWeightsForThisRGC;

        % Set the theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
        % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices

        centerConeIndicesAllRGCs{iSortedRGCindex} = centerConeIndicesForThisRGC(:);
        centerConeWeightsAllRGCs{iSortedRGCindex} = centerConeWeightsForThisRGC(:);
        centerRGCindicesAllRGCs{iSortedRGCindex} = repmat(iRGC, [numel(centerConeIndicesForThisRGC) 1]);

        surroundConeIndicesAllRGCs{iSortedRGCindex} = surroundConeIndicesForThisRGC(:);
        surroundConeWeightsAllRGCs{iSortedRGCindex} = surroundConeWeightsForThisRGC(:);
        surroundRGCindicesAllRGCs{iSortedRGCindex} = repmat(iRGC, [numel(surroundConeWeightsForThisRGC) 1]);


%         obj.rgcRFcenterConePoolingMatrix(pooledConeIndicesAndWeightsForThisRGC.centerConeIndices,iRGC) = ...
%             pooledConeIndicesAndWeightsForThisRGC.centerConeWeights;
%         obj.rgcRFsurroundConePoolingMatrix(pooledConeIndicesAndWeightsForThisRGC.surroundConeIndices,iRGC) = ...
%             pooledConeIndicesAndWeightsForThisRGC.surroundConeWeights;
%       
        if (abs((max(pooledConeIndicesAndWeightsForThisRGC.centerConeWeights(:)) - 1)) > 0.001)
            error('Max center weight must be 1. It is %2.4f', max(pooledConeIndicesAndWeightsForThisRGC.centerConeWeights(:)));
        end

    end % parfor iSortedRGCindex 

    centerConeIndices = vertcat(centerConeIndicesAllRGCs{:});
    centerConeWeights = vertcat(centerConeWeightsAllRGCs{:});
    rgcIndices = vertcat(centerRGCindicesAllRGCs{:});
    obj.rgcRFcenterConePoolingMatrix = sparse(centerConeIndices, rgcIndices, centerConeWeights);


    surroundConeIndices = vertcat(surroundConeIndicesAllRGCs{:});
    surroundConeWeights = vertcat(surroundConeWeightsAllRGCs{:});
    rgcIndices = vertcat(surroundRGCindicesAllRGCs{:});
    obj.rgcRFsurroundConePoolingMatrix = sparse(surroundConeIndices, rgcIndices, surroundConeWeights);

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
    totalCenterStrengthForModelRF, totalSurroundStrengthForModelRF)

    % Extract the retinal cone pooling weights compute function
    theRetinalConePoolingWeightsComputeFunction = theRFcomputeStruct.modelConstants.weightsComputeFunctionHandle;
    
    % Extract the retinal cone pooling weights params vector
    theRetinalConePoolingParamsVector = theRFcomputeStruct.retinalConePoolingParams.finalValues;

    % Extract the model constants
    modelConstants = theRFcomputeStruct.modelConstants;
       

    % Check to see if we compensate
    compensateForVariationsInConeEfficiency = ~(isnan(totalCenterStrengthForModelRF));

    if (compensateForVariationsInConeEfficiency)
        % Extract the fitted model center and surround retinal cone pooling weights
        fittedPooledConeIndicesAndWeights = theRetinalConePoolingWeightsComputeFunction(...
                modelConstants, theRetinalConePoolingParamsVector);
    
        fittedCenterConeWeights = fittedPooledConeIndicesAndWeights.centerConeWeights;
        fittedSurroundConeWeights = fittedPooledConeIndicesAndWeights.surroundConeWeights;
        fittedNonConnectableSurroundConeWeights = fittedPooledConeIndicesAndWeights.nonConnectableSurroundConeWeights;
    end


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