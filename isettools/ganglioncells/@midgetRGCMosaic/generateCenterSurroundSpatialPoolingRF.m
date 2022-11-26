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

    
    doCenterSurroundCorrections = false;
    if (doCenterSurroundCorrections)
        centerCorrectionFactors = zeros(1, rgcsNum);
        surroundCorrectionFactors = zeros(1, rgcsNum);
        % Initialize the total center/surround strengths for the model RFs
        totalCenterStrengthForModelRF = nan(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
        totalSurroundStrengthForModelRF = nan(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
    end

    % Sort RGCs with respect to their distance from the mRGC mosaic center
    distancesToEccGrid = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mean(obj.rgcRFpositionsDegs,1))).^2,2);
    [~,sortedRGCindices] = sort(distancesToEccGrid, 'ascend');


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

        fprintf('Computing c/s cone pooling weights for RGC %d of %d using the RTVFTobj #%d\n', iSortedRGCindex, rgcsNum, iObj);

        % Extract the retinal cone pooling weights compute function
        theRetinalConePoolingWeightsComputeFunction = theRTVFTobj.rfComputeStruct.modelConstants.weightsComputeFunctionHandle;
    
        % Extract the retinal cone pooling weights params vector
        theRetinalConePoolingParamsVector = theRTVFTobj.rfComputeStruct.retinalConePoolingParams.finalValues;

        % Extract the model constants
        modelConstants = theRTVFTobj.rfComputeStruct.modelConstants;
        
        % Update the center cone indices and weights for the target midgetRGC center 
        % so we can pass that info to theRetinalConePoolingWeightsComputeFunction()
        modelConstants.indicesOfCenterCones = indicesOfCenterCones;
        modelConstants.weightsOfCenterCones = connectivityVector(indicesOfCenterCones);
    
        % Compute the pooled cone indices and weights using the current
        % theRetinalConePoolingWeightsComputeFunction &
        % theRetinalConePoolingParamsVector
        pooledConeIndicesAndWeights = theRetinalConePoolingWeightsComputeFunction(...
            modelConstants, theRetinalConePoolingParamsVector);
    
        if (doCenterSurroundCorrections)
            % Here: compare pooledConeIndicesAndWeights (for this cell) to
            % the pooledConeIndicesAndWeights used to fit the RTV (i.e.
            % theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights) to
            % compute center and surround strength correction factors
            % what adjustments are needed. 

            if (isnan(totalCenterStrengthForModelRF(iObj)))
                idx = find(theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeWeights > 0);
                centerConeIndices = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeIndices(idx);
                relativeOuterSegmentLengths = obj.inputConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(centerConeIndices);
                relativeAbsorptionEfficacies = (obj.inputConeMosaic.coneApertureDiametersDegs(centerConeIndices)).^2 .* relativeOuterSegmentLengths;
                weights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.centerConeWeights(idx);
                
                totalCenterStrengthForModelRF(iObj) = sum(weights(:) .* relativeAbsorptionEfficacies(:),1)
              
                % Compute theRetinalRFcenterConeMap
                theConeApertureDiameterDegs = modelConstants.theConeMosaic.coneApertureDiametersDegs(centerConeIndices);
                theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(centerConeIndices,:);
                rfCenterPositionDegs = mean(theConePositionsDegs,1);
                theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
                theRetinalRFcenterConeMap{iObj} = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
                        modelConstants.theConeMosaic, ...
                        theConeApertureDiameterDegs, ...
                        theConePositionsDegs, ...
                        weights, ...
                        modelConstants.spatialSupportDegs);
            end

            if (isnan(totalSurroundStrengthForModelRF(iObj)))
                idx = find(theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeWeights > 0);
                surroundConeIndices = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeIndices(idx);
                relativeOuterSegmentLengths = obj.inputConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(surroundConeIndices);
                relativeAbsorptionEfficacies = (obj.inputConeMosaic.coneApertureDiametersDegs(surroundConeIndices)).^2 .* relativeOuterSegmentLengths;                                         
                weights = theRTVFTobj.rfComputeStruct.pooledConeIndicesAndWeights.surroundConeWeights(idx);
                totalSurroundStrengthForModelRF(iObj) = sum(weights(:) .* relativeAbsorptionEfficacies(:),1);
                surroundConesNumModel(iObj) = numel(idx);
    
                % Compute theRetinalRFsurroundConeMap
                theConeApertureDiameterDegs = modelConstants.theConeMosaic.coneApertureDiametersDegs(surroundConeIndices);                         
                theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(surroundConeIndices,:);
                theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
                theRetinalRFsurroundConeMap{iObj} = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
                    modelConstants.theConeMosaic, ...
                    theConeApertureDiameterDegs, ...
                    theConePositionsDegs, ...
                    weights, ...
                    modelConstants.spatialSupportDegs);
            end

            % Center strength correction factor
            idx = find(pooledConeIndicesAndWeights.centerConeWeights > 0);
            centerConeIndices = pooledConeIndicesAndWeights.centerConeIndices(idx)
            relativeOuterSegmentLengths = obj.inputConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(centerConeIndices);
            relativeAbsorptionEfficacies = (obj.inputConeMosaic.coneApertureDiametersDegs(centerConeIndices)).^2 .* relativeOuterSegmentLengths                   
            weights = pooledConeIndicesAndWeights.centerConeWeights(idx);
            totalCenterStrengthForThisRF = sum(weights(:) .* relativeAbsorptionEfficacies(:),1)
            centerStrengthCorrectionFactor = totalCenterStrengthForModelRF(iObj)/totalCenterStrengthForThisRF
            
            % Surround strength correction factor
            idx = find(pooledConeIndicesAndWeights.surroundConeWeights > 0);
            surroundConeIndices = pooledConeIndicesAndWeights.surroundConeIndices(idx);
            relativeOuterSegmentLengths = obj.inputConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(surroundConeIndices);
            relativeAbsorptionEfficacies = (obj.inputConeMosaic.coneApertureDiametersDegs(surroundConeIndices)).^2 .* relativeOuterSegmentLengths;                               
            weights = pooledConeIndicesAndWeights.surroundConeWeights(idx);
            totalSurroundStrengthForThisRF = sum(weights(:) .* relativeAbsorptionEfficacies(:),1);
            surroundStrengthCorrectionFactor = totalSurroundStrengthForModelRF(iObj)/totalSurroundStrengthForThisRF
        
     
            % Update theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
            % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices
            pooledConeIndicesAndWeights.centerConeWeights= ...
                pooledConeIndicesAndWeights.centerConeWeights * centerStrengthCorrectionFactor;

           pooledConeIndicesAndWeights.surroundConeWeights = ...
                pooledConeIndicesAndWeights.surroundConeWeights * surroundStrengthCorrectionFactor;

            % Keep track for printing the stats
            centerCorrectionFactors(iRGC) = (centerStrengthCorrectionFactor-1)*100;
            surroundCorrectionFactors(iRGC) = (surroundStrengthCorrectionFactor-1)*100;
        end % doCenterSurroundCorrections


        % Set the theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
        % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices
        obj.rgcRFcenterConePoolingMatrix(pooledConeIndicesAndWeights.centerConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.centerConeWeights;
        obj.rgcRFsurroundConePoolingMatrix(pooledConeIndicesAndWeights.surroundConeIndices,iRGC) = ...
            pooledConeIndicesAndWeights.surroundConeWeights;
       

        % Plot how the sum(surroundWeights)/sum(centerWeights) changes with
        % eccentricity within the patch
        figure(1001);
        if (iSortedRGCindex == 1)
            clf;
        end

        ax = subplot(1, numel(theRetinaToVisualFieldTransformerOBJList), iObj);
        hold(ax, 'on');

        plot(ax, distancesToEccGrids(iObj), sum(pooledConeIndicesAndWeights.surroundConeWeights) / sum(pooledConeIndicesAndWeights.centerConeWeights), 'k.');
        title(sprintf('RGC %d of %d', iSortedRGCindex,numel(sortedRGCindices)))
        
        xlabel(ax,'distance to nearest ecc grid (degs)');
        set(ax, 'XLim', [0 1], 'YLim', [0 1]);
        ylabel(ax, 'integrated S/C weights ratio')
        drawnow;
    end


    if (doCenterSurroundCorrections)
        fprintf('Center correction factors:   [ %+02.1f%% (min)  ..  %+02.1f%% (max)] \n', min(centerCorrectionFactors), max(centerCorrectionFactors));
        fprintf('Surround correction factors: [ %+02.1f%% (min)  ..  %+02.1f%% (max)] \n', min(surroundCorrectionFactors), max(surroundCorrectionFactors)); 
    end

  
    
end