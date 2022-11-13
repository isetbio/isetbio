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

    % Initialize the total center/surround strengths for the model RFs
    totalCenterStrengthForModelRF = nan(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));
    totalSurroundStrengthForModelRF = nan(1, numel(obj.theRetinaToVisualFieldTransformerOBJList));

    centerCorrectionFactors = zeros(1, rgcsNum);
    surroundCorrectionFactors = zeros(1, rgcsNum);

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

        fprintf('Computing c/s cone pooling weights for RGC %d of %d using the RTVFTobj #%d\n', iRGC, rgcsNum, iObj);
        
        theRTVFTobj
        theRTVFTobj.rfComputeStruct
        theRTVFTobj.rfComputeStruct.modelConstants
        pause
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
            totalCenterStrengthForModelRF(iObj) = sum(weights(:) .* relativeAbsorptionEfficacies(:),1);
            
            % Compute theRetinalRFcenterConeMap
            theConeAperturesDegs = modelConstants.theConeMosaic.coneApertureDiametersDegs(centerConeIndices);
    theConeCharacteristicRadiiDegs = ...
        modelConstants.coneCharacteristicRadiusConversionFactor * theConeAperturesDegs;
                                     
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(centerConeIndices,:);
    rfCenterPositionDegs = mean(theConePositionsDegs,1);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
    theRetinalRFcenterConeMap{iObj} = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeCharacteristicRadiiDegs, theConeAperturesDegs, theConePositionsDegs, weights, ...
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
             theConeAperturesDegs = modelConstants.theConeMosaic.coneApertureDiametersDegs(surroundConeIndices);
    theConeCharacteristicRadiiDegs = modelConstants.coneCharacteristicRadiusConversionFactor * theConeAperturesDegs;
                                     
    theConePositionsDegs = modelConstants.theConeMosaic.coneRFpositionsDegs(surroundConeIndices,:);
    theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
    theRetinalRFsurroundConeMap{iObj} = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
        theConeCharacteristicRadiiDegs, theConeAperturesDegs, theConePositionsDegs, weights, ...
        modelConstants.spatialSupportDegs);

        end

        % Center strength correction factor
        idx = find(pooledConeIndicesAndWeights.centerConeWeights > 0);
        centerConeIndices = pooledConeIndicesAndWeights.centerConeIndices(idx);
        relativeOuterSegmentLengths = obj.inputConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(centerConeIndices);
        relativeAbsorptionEfficacies = (obj.inputConeMosaic.coneApertureDiametersDegs(centerConeIndices)).^2 .* relativeOuterSegmentLengths;                      
        weights = pooledConeIndicesAndWeights.centerConeWeights(idx);
        totalCenterStrengthForThisRF = sum(weights(:) .* relativeAbsorptionEfficacies(:),1);
        centerStrengthCorrectionFactor = totalCenterStrengthForModelRF(iObj)/totalCenterStrengthForThisRF;

        % Surround strength correction factor
        idx = find(pooledConeIndicesAndWeights.surroundConeWeights > 0);
        surroundConeIndices = pooledConeIndicesAndWeights.surroundConeIndices(idx);
        relativeOuterSegmentLengths = obj.inputConeMosaic.outerSegmentLengthEccVariationAttenuationFactors(surroundConeIndices);
        relativeAbsorptionEfficacies = (obj.inputConeMosaic.coneApertureDiametersDegs(surroundConeIndices)).^2 .* relativeOuterSegmentLengths;                               
        weights = pooledConeIndicesAndWeights.surroundConeWeights(idx);
        totalSurroundStrengthForThisRF = sum(weights(:) .* relativeAbsorptionEfficacies(:),1);
        surroundStrengthCorrectionFactor = totalSurroundStrengthForModelRF(iObj)/totalSurroundStrengthForThisRF;

        % Update theMidgetRGCmosaic.rgcRFcenterConePoolingMatrix and
        % theMidgetRGCmosaic.rgcRFsurroundConePoolingMatrix sparse matrices
        centerIndices = pooledConeIndicesAndWeights.centerConeIndices;
        centerWeights = pooledConeIndicesAndWeights.centerConeWeights * centerStrengthCorrectionFactor;
        surroundIndices = pooledConeIndicesAndWeights.surroundConeIndices;
        surroundWeights = pooledConeIndicesAndWeights.surroundConeWeights * surroundStrengthCorrectionFactor;
        obj.rgcRFcenterConePoolingMatrix(centerIndices,iRGC) = centerWeights;
        obj.rgcRFsurroundConePoolingMatrix(surroundIndices,iRGC) = surroundWeights;

        actualModel = [sum(surroundWeights)/sum(centerWeights) sum(totalSurroundStrengthForModelRF(iObj))/totalCenterStrengthForModelRF(iObj)]
        contModel = [sum(sum(theRetinalRFsurroundConeMap{iObj}))/sum(sum(theRetinalRFcenterConeMap{iObj}))]
        [numel(surroundIndices) surroundConesNumModel(iObj) ]
    
        % Keep track for printing the stats
        centerCorrectionFactors(iRGC) = (centerStrengthCorrectionFactor-1)*100;
        surroundCorrectionFactors(iRGC) = (surroundStrengthCorrectionFactor-1)*100;
    end

    fprintf('Center correction factors:   [ %+02.1f%% (min)  ..  %+02.1f%% (max)] \n', min(centerCorrectionFactors), max(centerCorrectionFactors));
    fprintf('Surround correction factors: [ %+02.1f%% (min)  ..  %+02.1f%% (max)] \n', min(surroundCorrectionFactors), max(surroundCorrectionFactors)); 

end