function generateCenterSurroundSpatialPoolingRF(obj, theRetinaToVisualFieldTransformerOBJList, ...
            theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, ...
            theVisualSTFSurroundToCenterRcRatioGrid, ...
            theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid)

    obj.theRetinaToVisualFieldTransformerOBJList = theRetinaToVisualFieldTransformerOBJList;
    obj.theOpticsPositionGrid = theOpticsPositionGrid;
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
        centerCorrectionFactors = zeros(1, rgcsNum);
        surroundCorrectionFactors = zeros(1, rgcsNum);
        
        % Initialize the total center/surround strengths for the model RFs
        totalCenterStrengthForModelRF = nan(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
        totalSurroundStrengthForModelRF = nan(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
    end

    % Sort RGCs with respect to their distance from the mRGC mosaic center
    distancesToEccGrid = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mean(obj.rgcRFpositionsDegs,1))).^2,2);
    [~,sortedRGCindices] = sort(distancesToEccGrid, 'ascend');


    hFigModels = figure(1000); clf;
    set(hFigModels, 'Color', [1 1 1]);

    for iSortedRGCindex = 1:numel(sortedRGCindices)
    
        % Retrieve RGC index
        iRGC = sortedRGCindices(iSortedRGCindex);

        % Retrieve this cell's center cone indices and weights
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        indicesOfCenterCones = find(connectivityVector > 0.0001);

        % Compute distances of this RGC to all eccGrids
        distancesToEccGrids = sqrt(sum((bsxfun(@minus, obj.theOpticsPositionGrid, obj.rgcRFpositionsDegs(iRGC,:))).^2,2));

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
      

        % Visualize model fits
        if (isnan(totalCenterStrengthForModelRF(iObj)))
            figure(hFigModels);
            fprintf('Visualizing retinal maps for obj %d of %d\n', iObj, numel(totalCenterStrengthForModelRF));
            visualizedSpatialExtent = 0.1*[-1 1];
            ax = subplot(numel(obj.theRetinaToVisualFieldTransformerOBJList), 2, (iObj-1)*2+1);
            imagesc(ax, modelConstants.spatialSupportDegs(:,1), ...
                        modelConstants.spatialSupportDegs(:,2), ...
                        theRTVFTobj.rfComputeStruct.theRetinalRFcenterConeMap);
            axis(ax, 'image');
            set(ax, 'XLim', visualizedSpatialExtent, 'YLim', visualizedSpatialExtent, 'FontSize', 16);
    
            ax = subplot(numel(obj.theRetinaToVisualFieldTransformerOBJList), 2, (iObj-1)*2+2);
            imagesc(ax, modelConstants.spatialSupportDegs(:,1), ...
                        modelConstants.spatialSupportDegs(:,2), ...
                        theRTVFTobj.rfComputeStruct.theRetinalRFsurroundConeMap);
            axis(ax, 'image');
            set(ax, 'XLim', visualizedSpatialExtent, 'YLim', visualizedSpatialExtent, 'FontSize', 16);
    
            drawnow;
            colormap(gray(1024));
        end


        % Target and achieved ratios
        targetRsRcRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterRcRatio;
        targetSCintSensRatio = theRTVFTobj.targetVisualRFDoGparams.surroundToCenterIntegratedSensitivityRatio;
        fittedRsRcRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedRsRcRatio;
        fittedSCintSensRatio = theRTVFTobj.rfComputeStruct.theSTF.fittedSCIntSensRatio;

        
        fprintf('\n\n*** Stats for RGC %d ***\n', iSortedRGCindex);
        fprintf('Target Rs/Rc ratio: %2.2f, achieved: %2.2f\n', targetRsRcRatio, fittedRsRcRatio);
        fprintf('Target S/C int. sens. ratio: %2.2f, achieved: %2.2f\n', targetSCintSensRatio, fittedSCintSensRatio);
        fprintf('Computing c/s cone pooling weights for RGC %d of %d using the RTVFTobj #%d\n', iSortedRGCindex, rgcsNum, iObj);
        
        
        % Update the center cone indices and weights for the target midgetRGC center 
        % so we can pass that info to theRetinalConePoolingWeightsComputeFunction()
        modelConstants.indicesOfCenterCones = indicesOfCenterCones;
        modelConstants.weightsOfCenterCones = connectivityVector(indicesOfCenterCones);


        % Compute the pooled cone indices and weights using the
        % theRetinalConePoolingWeightsComputeFunction & theRetinalConePoolingParamsVector
        pooledConeIndicesAndWeights = theRetinalConePoolingWeightsComputeFunction(...
            modelConstants, theRetinalConePoolingParamsVector);


        % Extract the model center and surround retinal cone pooling weights
        currentRGCCenterConeWeights = pooledConeIndicesAndWeights.centerConeWeights;

        currentRGCSurroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights;
        currentRGCNonConnectableSurroundConeWeights = pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights;

        fprintf('fitted model RGC center/surround/surround non-connectable total weights: %2.2f (%d cones)/%2.2f (%d cones)/%2.2f (%d cones)\n', ...
            sum(fittedCenterConeWeights), numel(fittedCenterConeWeights), ...
            sum(fittedSurroundConeWeights), numel(fittedSurroundConeWeights), ...
            sum(fittedNonConnectableSurroundConeWeights), numel(fittedNonConnectableSurroundConeWeights));

        fprintf('current RGC center/surround total weights: %2.2f (%d cones)/%2.2f (%d cones)/%2.2f (%d cones)\n', ...
            sum(currentRGCCenterConeWeights), numel(currentRGCCenterConeWeights), ...
            sum(currentRGCSurroundConeWeights), numel(currentRGCSurroundConeWeights), ...
            sum(currentRGCNonConnectableSurroundConeWeights), numel(currentRGCNonConnectableSurroundConeWeights));
        
        if (compensateForVariationsInConeEfficiency)
            % Here: compare pooledConeIndicesAndWeights (for this cell) to
            % the pooledConeIndicesAndWeights used to fit the RTV (i.e.
            % theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights) to
            % compute center and surround strength correction factors
            % what adjustments are needed. 

            if (isnan(totalCenterStrengthForModelRF(iObj)))
                idx = find(theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeWeights > 0);
                weights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeWeights(idx);
                totalCenterStrengthForModelRF(iObj) = sum(weights(:));
            end

            if (isnan(totalSurroundStrengthForModelRF(iObj)))
                idx = find(theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeWeights > 0);
                weights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeWeights(idx);
                totalSurroundStrengthForModelRF(iObj) = sum(weights(:));
            end

            % Center cone indices and weights for the current RGC
            totalCenterStrengthForThisRF = sum(pooledConeIndicesAndWeights.centerConeWeights);
            centerStrengthCorrectionFactor = totalCenterStrengthForModelRF(iObj)/totalCenterStrengthForThisRF;
            pooledConeIndicesAndWeights.centerConeWeights = pooledConeIndicesAndWeights.centerConeWeights * centerStrengthCorrectionFactor;


%             % Surround strength correction factor
%             % First adjust for different # of S-cones
            if (1==1)
                fittedBoost = 1 + sum(fittedNonConnectableSurroundConeWeights)/sum(fittedSurroundConeWeights);
                currentBoost = 1 + sum(currentRGCNonConnectableSurroundConeWeights)/sum(currentRGCSurroundConeWeights);
                surroundBoostRatio = currentBoost/fittedBoost
                pooledConeIndicesAndWeights.surroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights * surroundBoostRatio;
            end


            totalSurroundStrengthForThisRF = sum(pooledConeIndicesAndWeights.surroundConeWeights);
            surroundStrengthCorrectionFactor = totalSurroundStrengthForModelRF(iObj)/totalSurroundStrengthForThisRF;
            pooledConeIndicesAndWeights.surroundConeWeights = pooledConeIndicesAndWeights.surroundConeWeights * surroundStrengthCorrectionFactor;


            % Keep track for printing the stats
            centerCorrectionFactors(iRGC) = (centerStrengthCorrectionFactor-1)*100;
            surroundCorrectionFactors(iRGC) = (surroundStrengthCorrectionFactor-1)*100;
        end % compensateForVariationsInConeEfficiency


        % Set the theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
        % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices
        obj.rgcRFcenterConePoolingMatrix(pooledConeIndicesAndWeights.centerConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.centerConeWeights;
        obj.rgcRFsurroundConePoolingMatrix(pooledConeIndicesAndWeights.surroundConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.surroundConeWeights;
       

        % Plot how the sum(surroundWeights)/sum(centerWeights) changes with
        % eccentricity within the patch
        hFig = figure(1001);
        if (iSortedRGCindex == 1)
            set(hFig, 'Color', [1 1 1]);
            clf;
        end

        ax = subplot(1, numel(theRetinaToVisualFieldTransformerOBJList), iObj);
        hold(ax, 'on');

        plot(ax, distancesToEccGrids(iObj), sum(pooledConeIndicesAndWeights.surroundConeWeights) / sum(pooledConeIndicesAndWeights.centerConeWeights), 'k.');
        title(sprintf('RGC %d of %d (%d center cones)', iSortedRGCindex,numel(sortedRGCindices), obj.theConesNumPooledByTheRFcenterGrid(iObj)))
        
        xlabel(ax,'distance to nearest ecc grid (degs)');
        set(ax, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
        ylabel(ax, 'integrated S/C weights ratio (retinal)')
        drawnow;
    end


    if (compensateForVariationsInConeEfficiency)
        fprintf('Center correction factors:   [ %+02.1f%% (min)  ..  %+02.1f%% (max)] \n', min(centerCorrectionFactors), max(centerCorrectionFactors));
        fprintf('Surround correction factors: [ %+02.1f%% (min)  ..  %+02.1f%% (max)] \n', min(surroundCorrectionFactors), max(surroundCorrectionFactors)); 
    end

  
    
end