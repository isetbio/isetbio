function generateComputeReadyMidgetRGCMosaic(obj, optimizedRGCpoolingObjectsFileName, ...
    computeReadyMosaicFilename, visualizeInterpolation)

    % Load the optimizedRFcomputeStructs (retinal cone pooling models)
    % across all grid nodes in the mosaic
    [LconeRFcomputeStructsList, MconeRFcomputeStructsList] = ...
        loadOptimizedRFcomputeStructs(obj, optimizedRGCpoolingObjectsFileName);

    % Compute indices of cones connected to the centers and the surrounds
    % and the weights of these connections for all cells in theRGCMosaic
    [centerConeIndicesAllRGCs, centerConeWeightsAllRGCs, centerRGCindicesAllRGCs, ...
     surroundConeIndicesAllRGCs, surroundConeWeightsAllRGCs, surroundRGCindicesAllRGCs ] = ...
         computePoolingWeightsByInterpolatingModelsAcrossTheSamplingGrid(...
                obj, LconeRFcomputeStructsList, MconeRFcomputeStructsList, ...
                visualizeInterpolation);

    % Generate sparse matrix encoding cone connections to the RF centers
    centerConeIndices = vertcat(centerConeIndicesAllRGCs{:});
    centerConeWeights = vertcat(centerConeWeightsAllRGCs{:});
    rgcIndices = vertcat(centerRGCindicesAllRGCs{:});
    rgcRFcenterConePoolingMatrix = sparse(centerConeIndices, rgcIndices, centerConeWeights, obj.theRGCMosaic.inputConeMosaic.conesNum, obj.theRGCMosaic.rgcsNum);

    % Sanity check
    sanityCheck(rgcRFcenterConePoolingMatrix, 'center');

    % Generate sparse matrix encoding cone connections to the RF surrounds
    surroundConeIndices = vertcat(surroundConeIndicesAllRGCs{:});
    surroundConeWeights = vertcat(surroundConeWeightsAllRGCs{:});
    rgcIndices = vertcat(surroundRGCindicesAllRGCs{:});
    rgcRFsurroundConePoolingMatrix = sparse(surroundConeIndices, rgcIndices, surroundConeWeights, obj.theRGCMosaic.inputConeMosaic.conesNum, obj.theRGCMosaic.rgcsNum);

    % Sanity check
    sanityCheck(rgcRFsurroundConePoolingMatrix, 'surround');

    % All good. Proceed with baking in the computed center/surround cone pooling matrices to
    % theRCMosaic. This also sets the center-connectivity sparse matrix to []. 
    obj.theRGCMosaic.bakeInConePoolingMatrices(...
        rgcRFcenterConePoolingMatrix, rgcRFsurroundConePoolingMatrix);


    

    % Export this compute-ready midgetRGCMosaic
    theComputeReadyMRGCmosaic = obj.theRGCMosaic;
    save(computeReadyMosaicFilename, 'theComputeReadyMRGCmosaic', '-v7.3');
    fprintf('The compute-ready mRGCMosaic was exported in %s.\n', computeReadyMosaicFilename);
end

function [centerConeIndicesAllRGCs, centerConeWeightsAllRGCs, centerRGCindicesAllRGCs, ...
          surroundConeIndicesAllRGCs, surroundConeWeightsAllRGCs, surroundRGCindicesAllRGCs ] = ...
                computePoolingWeightsByInterpolatingModelsAcrossTheSamplingGrid(obj, ...
                    LconeRFcomputeStructsList, MconeRFcomputeStructsList, visualizeInterpolation)

    rgcsNum = size(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix,2);
   
    centerConeIndicesAllRGCs = cell(rgcsNum,1);
    centerConeWeightsAllRGCs = cell(rgcsNum,1);
    centerRGCindicesAllRGCs = cell(rgcsNum,1);

    surroundConeIndicesAllRGCs = cell(rgcsNum,1);
    surroundConeWeightsAllRGCs = cell(rgcsNum,1);
    surroundRGCindicesAllRGCs = cell(rgcsNum,1);

    if (visualizeInterpolation)
        for theCurrentRGCindex = 1:rgcsNum 
            fprintf('Generating surround for RGC%d of %d RGC\n', theCurrentRGCindex, rgcsNum);

            d = computeConeInputsForRCG(obj, theCurrentRGCindex, ...
                LconeRFcomputeStructsList, MconeRFcomputeStructsList, true);

            centerConeIndicesAllRGCs{theCurrentRGCindex} = d.centerConeIndices;
            centerConeWeightsAllRGCs{theCurrentRGCindex} = d.centerConeWeights;
            centerRGCindicesAllRGCs{theCurrentRGCindex} = repmat(theCurrentRGCindex, [numel(d.centerConeIndices) 1]);

            surroundConeIndicesAllRGCs{theCurrentRGCindex} = d.surroundConeIndices;
            surroundConeWeightsAllRGCs{theCurrentRGCindex} = d.surroundConeWeights;
            surroundRGCindicesAllRGCs{theCurrentRGCindex} = repmat(theCurrentRGCindex, [numel(d.surroundConeIndices) 1]);
        end
    else
        parfor theCurrentRGCindex = 1:rgcsNum 
            fprintf('Generating surround for RGC%d of %d RGC\n', theCurrentRGCindex, rgcsNum);

            d = computeConeInputsForRCG(obj, theCurrentRGCindex, ...
                LconeRFcomputeStructsList, MconeRFcomputeStructsList, false);

            centerConeIndicesAllRGCs{theCurrentRGCindex} = d.centerConeIndices;
            centerConeWeightsAllRGCs{theCurrentRGCindex} = d.centerConeWeights;
            centerRGCindicesAllRGCs{theCurrentRGCindex} = repmat(theCurrentRGCindex, [numel(d.centerConeIndices) 1]);

            surroundConeIndicesAllRGCs{theCurrentRGCindex} = d.surroundConeIndices;
            surroundConeWeightsAllRGCs{theCurrentRGCindex} = d.surroundConeWeights;
            surroundRGCindicesAllRGCs{theCurrentRGCindex} = repmat(theCurrentRGCindex, [numel(d.surroundConeIndices) 1]);
        end
    end
end



function dStruct = computeConeInputsForRCG(obj, theCurrentRGCindex, ...
    LconeRFcomputeStructsList, MconeRFcomputeStructsList, visualizeInterpolation)

    % The cell position
    theCurrentRGCposition = obj.theRGCMosaic.rgcRFpositionsDegs(theCurrentRGCindex,:);

    % Retrieve this cell's # of center cone indices
    connectivityVector = full(squeeze(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:, theCurrentRGCindex)));
    centerConeIndices = find(connectivityVector > 0.0001);

    % Retrieve this cell's center cone weights
    weightsOfCenterConesForCurrentRGC = connectivityVector(centerConeIndices);

    % Determine the net L-cone, and net M-cone weights in the RF center of this RGC
    [theCenterConeTypeWeights, ~, ~, theCenterConeTypes] = obj.centerConeTypeWeights(theCurrentRGCindex);
    theTotalCenterConeWeights = sum(theCenterConeTypeWeights(:));

    % Spectral interpolation: compute the weights for L- and M-cone compute structs to be proportional to the
    % # of center L- and M-cones, respectively
    theLconeComputeStructWeight = theCenterConeTypeWeights(find(theCenterConeTypes == cMosaic.LCONE_ID)) / theTotalCenterConeWeights;
    theMconeComputeStructWeight = theCenterConeTypeWeights(find(theCenterConeTypes == cMosaic.MCONE_ID)) / theTotalCenterConeWeights;

    % Find all model RGCs with the same number of center cones as this RGC
    gridNodeIndicesWithMatchedCenterConesNum = find(obj.conesNumPooledByTheRFcenterGrid == numel(centerConeIndices));

    % Find the model RGC indices with this many L-cones in their RF center
    targetLconeModelRGCindices = obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndicesWithMatchedCenterConesNum);
    targetLconeModelRGCpositions = obj.theRGCMosaic.rgcRFpositionsDegs(targetLconeModelRGCindices,:);
    
    % Spatial interpolation for L-center models: compute this RGC' weights for pooling parameters from the 3 triangularing L-center RGC models
    [triangulatingLconeModelRGCIndices, triangulatingLconeModelRGCweights, triangulatingLconeGridIndices] = ...
            obj.triangulatingGridNodeIndicesAndWeights(theCurrentRGCposition, ...
                targetLconeModelRGCpositions, targetLconeModelRGCindices, gridNodeIndicesWithMatchedCenterConesNum);

    % Find the model RGC indices with this many M-cones in their RF center
    targetMconeModelRGCindices = obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndicesWithMatchedCenterConesNum);
    targetMconeModelRGCpositions = obj.theRGCMosaic.rgcRFpositionsDegs(targetMconeModelRGCindices,:);

    % Spatial interpolation for M-center models: compute this RGC' weights for pooling parameters from the 3 triangularing M-center RGC models
    [triangulatingMconeModelRGCIndices, triangulatingMconeModelRGCweights, triangulatingMconeGridIndices] = ...
            obj.triangulatingGridNodeIndicesAndWeights(theCurrentRGCposition, ...
                targetMconeModelRGCpositions, targetMconeModelRGCindices, gridNodeIndicesWithMatchedCenterConesNum);

    % Spatio-spectral intepolation: combine spatial and spectral weights.
    triangulatingLconeModelRGCweights = triangulatingLconeModelRGCweights * theLconeComputeStructWeight;
    triangulatingMconeModelRGCweights = triangulatingMconeModelRGCweights * theMconeComputeStructWeight;

    % Update modelConstants for the current RGC
    modelConstants = LconeRFcomputeStructsList{triangulatingLconeGridIndices(1)}.modelConstants;
    modelConstantsForCurrentRGC = updateModelConstantsForCurrentRGC(obj.theRGCMosaic, modelConstants, ...
                theCurrentRGCposition, centerConeIndices, weightsOfCenterConesForCurrentRGC);

    coneIndicesAndWeightsForLconeModels = coneIndicesAndWeightsFromTriangulatingModels(...
            modelConstantsForCurrentRGC,...
            LconeRFcomputeStructsList, ...
            triangulatingLconeGridIndices, ...
            triangulatingLconeModelRGCweights);

       
    coneIndicesAndWeightsForMconeModels = coneIndicesAndWeightsFromTriangulatingModels(...
            modelConstantsForCurrentRGC,...
            MconeRFcomputeStructsList, ...
            triangulatingMconeGridIndices, ...
            triangulatingMconeModelRGCweights);

    % Accumulate center cone weights from all models now
    centerConeWeights = zeros(size(centerConeIndices));
    surroundConeIndices = [];
    surroundConeWeights = [];

    for iModel = 1:numel(coneIndicesAndWeightsForLconeModels)
        % Center cones
        centerConeWeights = centerConeWeights + ...
                    coneIndicesAndWeightsForLconeModels{iModel}.centerConeWeights;

        % Surround cones
        theModelSurroundConeWeights = coneIndicesAndWeightsForLconeModels{iModel}.surroundConeWeights;
        theModelSurroundConeIndices = coneIndicesAndWeightsForLconeModels{iModel}.surroundConeIndices;
        [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
                surroundConeIndices, surroundConeWeights, ...
                theModelSurroundConeIndices, theModelSurroundConeWeights);
    end

    for iModel = 1:numel(coneIndicesAndWeightsForMconeModels)
        % Center cones
        centerConeWeights = centerConeWeights + ...
                    coneIndicesAndWeightsForMconeModels{iModel}.centerConeWeights;

        % Surround cones
        theModelSurroundConeWeights = coneIndicesAndWeightsForMconeModels{iModel}.surroundConeWeights;
        theModelSurroundConeIndices = coneIndicesAndWeightsForMconeModels{iModel}.surroundConeIndices;
        [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
                surroundConeIndices, surroundConeWeights, ...
                theModelSurroundConeIndices, theModelSurroundConeWeights);

    end

    % Assemble struct with center and surround cone indices and weights
    dStruct = struct();
    dStruct.centerConeIndices = centerConeIndices(:);
    dStruct.centerConeWeights = centerConeWeights(:);
    dStruct.surroundConeIndices = surroundConeIndices(:);
    dStruct.surroundConeWeights = surroundConeWeights(:);

    if (visualizeInterpolation)
       MosaicPoolingOptimizer.visualizeMultiSpectralMultiFocalRFgeneration(1000, ...
                obj.theRGCMosaic.inputConeMosaic, theCurrentRGCindex, theCurrentRGCposition, ...
                cMosaic.LCONE_ID, coneIndicesAndWeightsForLconeModels, ...
                obj.theRGCMosaic.rgcRFpositionsDegs(obj.targetRGCindicesWithLconeMajorityCenter(triangulatingLconeGridIndices),:), ...
                triangulatingLconeModelRGCIndices, triangulatingLconeModelRGCweights, ...
                centerConeIndices, centerConeWeights, ...
                surroundConeIndices, surroundConeWeights);

       MosaicPoolingOptimizer.visualizeMultiSpectralMultiFocalRFgeneration(1000, ...
                obj.theRGCMosaic.inputConeMosaic, theCurrentRGCindex, theCurrentRGCposition, ...
                cMosaic.MCONE_ID, coneIndicesAndWeightsForMconeModels, ...
                obj.theRGCMosaic.rgcRFpositionsDegs(obj.targetRGCindicesWithMconeMajorityCenter(triangulatingMconeGridIndices),:), ...
                triangulatingMconeModelRGCIndices, triangulatingMconeModelRGCweights, ...
                centerConeIndices, centerConeWeights, ...
                surroundConeIndices, surroundConeWeights);

       pause
   end

end


function [surroundConeIndices, surroundConeWeights] = accumulateSurroundWeights(...
                surroundConeIndices, surroundConeWeights, ...
                newSurroundConeIndices, newSurroundConeWeights)

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


function pooledConeIndicesAndWeightsForModels = coneIndicesAndWeightsFromTriangulatingModels(...
    modelConstantsForCurrentRGC, theRFcomputeStructsList, ...
    triangulatingGridIndices, triangulatingGridWeights)

    pooledConeIndicesAndWeightsForModels = cell(1,numel(triangulatingGridIndices));

    for iModel = 1:numel(triangulatingGridIndices)
       % Retrieve theRFcomputeStruct for the model RGC
       theRFcomputeStruct = theRFcomputeStructsList{triangulatingGridIndices(iModel)};

       % Update theRFcomputeStruct with modelConstants for the current RGC
       theRFcomputeStruct.modelConstants = modelConstantsForCurrentRGC; 

       % Accumulate pooled cone indices and weights from this model
       % applied to the current RGC
       pooledConeIndicesAndWeightsForModels{iModel} = centerNormalizedPooledConeIndicesAndWeights(...
          theRFcomputeStruct, triangulatingGridWeights(iModel), ...
          theRFcomputeStruct.modelTotalCenterStrength, ...
          theRFcomputeStruct.modelTotalSurroundStrength);
     end % iModel
end


function sanityCheck(subregionConePoolingMatrix, subregionName)
    totalSubregionWeightsForAllRGCs= sum(subregionConePoolingMatrix,1);
    connectedRGCindices = numel(find(totalSubregionWeightsForAllRGCs > 0.0001));
    assert(connectedRGCindices == size(subregionConePoolingMatrix,2), ...
        sprintf('Some RGCs do not have inputs to their RF %s', subregionName));
end


function [LconeComputeStructsList, MconeComputeStructsList] = loadOptimizedRFcomputeStructs(obj, optimizedRGCpoolingObjectsFileName)

    LconeComputeStructsList = cell(1, numel(obj.conesNumPooledByTheRFcenterGrid));
    MconeComputeStructsList = cell(1, numel(obj.conesNumPooledByTheRFcenterGrid));

    for gridNodeIndex = 1:numel(obj.conesNumPooledByTheRFcenterGrid)
        % Optimized RGCpooling object filename
        optimizedRGCpoolingObjectsFileNameForThisNode = ...
            strrep(optimizedRGCpoolingObjectsFileName, '.mat', sprintf('_ForGridNode_%d.mat', gridNodeIndex));

        fprintf('Loading optimized L- and M-cone RF compute structs from node %d of %d (file: %s)\n', ...
            gridNodeIndex, numel(obj.conesNumPooledByTheRFcenterGrid), optimizedRGCpoolingObjectsFileNameForThisNode);
        load(optimizedRGCpoolingObjectsFileNameForThisNode, ...
            'theLconeRFcomputeStruct', ...
            'theMconeRFcomputeStruct');

        % Compute the L-cone models' total center/surround and non-connectable
        % surround strengths
        pooledConeIndicesAndWeights = centerNormalizedPooledConeIndicesAndWeights(...
                theLconeRFcomputeStruct, 1.0, [], []);
        theLconeRFcomputeStruct.modelTotalCenterStrength = sum(pooledConeIndicesAndWeights.centerConeWeights);
        theLconeRFcomputeStruct.modelTotalSurroundStrength = sum(pooledConeIndicesAndWeights.surroundConeWeights);
        theLconeRFcomputeStruct.modelTotalNonConnectableSurroundConesStrength = sum(pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights);


        % Compute the M-cone models' total center/surround and non-connectable
        % surround strengths
        pooledConeIndicesAndWeights = centerNormalizedPooledConeIndicesAndWeights(...
                theMconeRFcomputeStruct, 1.0, [], []);
        theMconeRFcomputeStruct.modelTotalCenterStrength = sum(pooledConeIndicesAndWeights.centerConeWeights);
        theMconeRFcomputeStruct.modelTotalSurroundStrength = sum(pooledConeIndicesAndWeights.surroundConeWeights);
        theMconeRFcomputeStruct.modelTotalNonConnectableSurroundConesStrength = sum(pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights);

        fprintf('L-cone model at grid node%d has total strengths: %2.3f (center) %2.3f (surround)\n', gridNodeIndex, ...
            theLconeRFcomputeStruct.modelTotalCenterStrength , theLconeRFcomputeStruct.modelTotalSurroundStrength );
       
        fprintf('M-cone model at grid node%d has total strengths: %2.3f (center) %2.3f (surround)\n', gridNodeIndex, ...
            theMconeRFcomputeStruct.modelTotalCenterStrength , theMconeRFcomputeStruct.modelTotalSurroundStrength );
        
        % Add to the list
        LconeComputeStructsList{gridNodeIndex} = theLconeRFcomputeStruct;
        MconeComputeStructsList{gridNodeIndex} = theMconeRFcomputeStruct;
    end % iNode

end

function updatedModelConstants = updateModelConstantsForCurrentRGC(...
    theRGCMosaic, modelConstants, theRGCpositionDegs, ...
    indicesOfConesPooledByTheRFcenter,weightsOfConesPooledByTheRFcenter)

    updatedModelConstants = modelConstants;
    
    % Maximum support for the surround, in degrees, taken as MosaicPoolingOptimizer.maxSurroundSupportFactor times the C&K
    % surrounds at the cell's eccentricity
    radialEccentricityForThisRGC = sqrt(sum(theRGCpositionDegs.^2,2));
    updatedModelConstants.maxSurroundSupportDegs = MosaicPoolingOptimizer.maxSurroundSupportFactor * ...
                RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(radialEccentricityForThisRGC);

    % Compute the RF center position
    RFcenterPos = mean(theRGCMosaic.inputConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);

    % Compute the distances of ALL cones in the input cone mosaic from the RF center
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, theRGCMosaic.inputConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);
    coneDistancesFromRFCenter = sqrt(coneDistancesFromRFCenterSquared);
    updatedModelConstants.cachedData.surroundConeIndices = find(coneDistancesFromRFCenter <= updatedModelConstants.maxSurroundSupportDegs);
    updatedModelConstants.cachedData.coneDistancesFromRFCenter = coneDistancesFromRFCenter(updatedModelConstants.cachedData.surroundConeIndices);

    updatedModelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
    updatedModelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;
end


function pooledConeIndicesAndWeights = centerNormalizedPooledConeIndicesAndWeights(...
    theRFcomputeStruct, triangulatingWeight, modelTotalCenterStrength, modelTotalSurroundStrength)

     thePoolingWeightsComputeFunction = theRFcomputeStruct.modelConstants.weightsComputeFunctionHandle;
     optimizedRetinalPoolingParams = theRFcomputeStruct.retinalConePoolingParams.finalValues;

     % Compute pooled cone indices and weights
     pooledConeIndicesAndWeights = thePoolingWeightsComputeFunction(...
         theRFcomputeStruct.modelConstants, optimizedRetinalPoolingParams);

     % Compute the max center cone weight
     maxCenterWeight = max(pooledConeIndicesAndWeights.centerConeWeights(:));

     % Divide all weights with the max center weight
     pooledConeIndicesAndWeights.centerConeWeights = ...
         pooledConeIndicesAndWeights.centerConeWeights / maxCenterWeight;
     pooledConeIndicesAndWeights.surroundConeWeights = ...
         pooledConeIndicesAndWeights.surroundConeWeights / maxCenterWeight;
     pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights = ...
         pooledConeIndicesAndWeights.nonConnectableSurroundConeWeights / maxCenterWeight;


     % Adjust for any differences in the total strength of the
     % center for the currentRGC vs the model
     if (~isempty(modelTotalCenterStrength))
         totalCenterStrengthForTheCurrentRGC = sum(pooledConeIndicesAndWeights.centerConeWeights);
         correctionFactor = modelTotalCenterStrength/totalCenterStrengthForTheCurrentRGC;
         pooledConeIndicesAndWeights.centerConeWeights = ...
             pooledConeIndicesAndWeights.centerConeWeights * correctionFactor;
     end

     % Adjust for any differences in the total strength of the
     % surround for the currentRGC vs the model
     if (~isempty(modelTotalSurroundStrength))
         totalSurroundStrengthForTheCurrentRGC = sum(pooledConeIndicesAndWeights.surroundConeWeights);
         correctionFactor = modelTotalSurroundStrength/totalSurroundStrengthForTheCurrentRGC;
         pooledConeIndicesAndWeights.surroundConeWeights = ...
             pooledConeIndicesAndWeights.surroundConeWeights * correctionFactor;
     end

     % Finally apply triangulating weight
     pooledConeIndicesAndWeights.centerConeWeights = ...
         triangulatingWeight * pooledConeIndicesAndWeights.centerConeWeights;
     pooledConeIndicesAndWeights.surroundConeWeights = ...
         triangulatingWeight * pooledConeIndicesAndWeights.surroundConeWeights;
end