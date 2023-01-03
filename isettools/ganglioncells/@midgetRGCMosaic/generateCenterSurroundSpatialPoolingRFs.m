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

    visualizedSpatialExtent = 0.1*[-1 1];

    if (1==2)
    hFigModels = figure(1000); clf;
    set(hFigModels, 'Color', [1 1 1]);

    
    fprintf('Visualizing retinal maps for obj %d of %d\n', iObj, numel(totalCenterStrengthForModelRF));
    
    for iObj = 1:numel(obj.theRetinaToVisualFieldTransformerOBJList)

        % Target and achieved ratios
        targetRsRcRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
        targetSCintSensRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
        fittedRsRcRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedRsRcRatio;
        fittedSCintSensRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedSCIntSensRatio;

        fprintf('Target Rs/Rc ratio: %2.2f, achieved: %2.2f\n', targetRsRcRatio, fittedRsRcRatio);
        fprintf('Target S/C int. sens. ratio: %2.2f, achieved: %2.2f\n', targetSCintSensRatio, fittedSCintSensRatio);

        theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};
        ax = subplot(2,numel(obj.theRetinaToVisualFieldTransformerOBJList), iObj);
        imagesc(ax, theRTVFTobj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1), ...
                    theRTVFTobj.rfComputeStruct.modelConstants.spatialSupportDegs(:,2), ...
                    theRTVFTobj.rfComputeStruct.theRetinalRFcenterConeMap);
        axis(ax, 'image');
        set(ax, 'XLim', visualizedSpatialExtent, 'YLim', visualizedSpatialExtent, 'FontSize', 16);

        ax = subplot(2,numel(obj.theRetinaToVisualFieldTransformerOBJList), numel(obj.theRetinaToVisualFieldTransformerOBJList)+iObj);
        imagesc(ax, theRTVFTobj.rfComputeStruct.modelConstants.spatialSupportDegs(:,1), ...
                    theRTVFTobj.rfComputeStruct.modelConstants.spatialSupportDegs(:,2), ...
                    theRTVFTobj.rfComputeStruct.theRetinalRFsurroundConeMap);
        axis(ax, 'image');
        set(ax, 'XLim', visualizedSpatialExtent, 'YLim', visualizedSpatialExtent, 'FontSize', 16);

        drawnow;
        colormap(gray(1024));
    end
    pause
    end


    % Sort RGCs with respect to their distance from the mRGC mosaic center
    distancesToEccGrid = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mean(obj.rgcRFpositionsDegs,1))).^2,2);
    [~,sortedRGCindices] = sort(distancesToEccGrid, 'ascend');


    for iSortedRGCindex = 1:numel(sortedRGCindices)
        fprintf('Generating surround for %d of %d RGC\n', iSortedRGCindex, numel(sortedRGCindices));

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
            fprintf('nearby RGC %d center weights: max = %f, sum = %f\n', ...
              max(pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.centerConeWeights), ...
              sum(pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.centerConeWeights) ...
              );
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

        visualizeRetinalRFsurroundMapFromWeigtedAcrossRTVFobjects = ~true;
        if (visualizeRetinalRFsurroundMapFromWeigtedAcrossRTVFobjects)
            figure(5678); clf;

            modelConstants = obj.theRetinaToVisualFieldTransformerOBJList{triangulatingRTVFobjIndices(1)}.rfComputeStruct.modelConstants;
            xx = modelConstants.spatialSupportDegs(:,1);
            idx = find(abs(xx)< max(visualizedSpatialExtent));
            spatialSupportDegs = modelConstants.spatialSupportDegs(idx,:);

            
            theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(pooledConeIndicesAndWeights.centerConeIndices,:);
            rfCenterPositionDegs = mean(theConePositionsDegs,1);
    
            for iNearbyObj = 1:numel(triangulatingRTVFobjIndices)
                ax = subplot(2,2,iNearbyObj);
                surroundConeIndices = pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.surroundConeIndices;
                surroundConeWeights = pooledConeIndicesAndWeightsForNearbyObj{iNearbyObj}.surroundConeWeights;
                theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(surroundConeIndices,:);
                theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
                theModelRetinalRFsurroundConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
                    modelConstants.theConeMosaic,...
                    theConePositionsDegs, ...
                    surroundConeWeights, ...
                    spatialSupportDegs);

                imagesc(ax,modelConstants.spatialSupportDegs(:,1), modelConstants.spatialSupportDegs(:,2), theModelRetinalRFsurroundConeMap);
                axis(ax, 'image');
                set(ax, 'XLim', visualizedSpatialExtent, 'YLim', visualizedSpatialExtent, 'FontSize', 16);
                title(ax,sprintf('nearby model %d surround (weight: %2.2f)\n(%d surround cones, %d center cones)', ...
                    iNearbyObj, triangulatingRTVFobjWeights(iNearbyObj), numel(surroundConeIndices), numel(pooledConeIndicesAndWeights.centerConeIndices)));
            end


            ax = subplot(2,2,4);
            theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(pooledConeIndicesAndWeights.surroundConeIndices,:);
            theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
            theRGCRetinalRFsurroundConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
                    modelConstants.theConeMosaic,...
                    theConePositionsDegs, ...
                    pooledConeIndicesAndWeights.surroundConeWeights, ...
                    spatialSupportDegs);

            imagesc(ax,modelConstants.spatialSupportDegs(:,1), modelConstants.spatialSupportDegs(:,2), theRGCRetinalRFsurroundConeMap);
            axis(ax, 'image')
            set(ax, 'XLim', visualizedSpatialExtent, 'YLim', visualizedSpatialExtent, 'FontSize', 16);
            title(ax, sprintf('RGC %d surround\n(%d surround cones)', iSortedRGCindex, numel(pooledConeIndicesAndWeights.surroundConeWeights)));
            colormap(gray(1024))
            drawnow;
        end



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

    beVerbose = false;
    if (beVerbose)
    fprintf('fitted model RGC center/surround/surround/non-connectable total weights: %2.2f (%d cones)/%2.2f (%d cones)/%2.2f (%d S-cones)\n', ...
            sum(fittedCenterConeWeights), numel(fittedCenterConeWeights), ...
            sum(fittedSurroundConeWeights), numel(fittedSurroundConeWeights), ...
            sum(fittedNonConnectableSurroundConeWeights), numel(fittedNonConnectableSurroundConeWeights));

    fprintf('current RGC center/surround/non-connectable total weights: %2.2f (%d cones)/%2.2f (%d cones)/%2.2f (%d S-cones)\n', ...
            sum(currentRGCCenterConeWeights), numel(currentRGCCenterConeWeights), ...
            sum(currentRGCSurroundConeWeights), numel(currentRGCSurroundConeWeights), ...
            sum(currentRGCNonConnectableSurroundConeWeights), numel(currentRGCNonConnectableSurroundConeWeights));
    end

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