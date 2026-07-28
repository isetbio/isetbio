function cachedData = cache(theMRGCMosaic, ...
  theTargetRGCindex, poolingOptimizationParamsStruct)

  % Indices of RF center input cones
  % Get all center cones, and tell the method to not warn us about difference in cone input numerosity
  minConeWeightIncluded = 0.001;
  inputConeIndices = theMRGCMosaic.singleCellConnectivityStats(...
    theTargetRGCindex, 'center', ...
    'minConeWeightIncluded', minConeWeightIncluded, ...
    'inputConeIndicesOnly', true, ...
    'warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', false);

  cachedData.centerConeWeights  = full(theMRGCMosaic.rgcRFcenterConeConnectivityMatrix(inputConeIndices, theTargetRGCindex));
  cachedData.centerConeIndices = inputConeIndices;

  % Compute the RF center position
  cachedData.RFcenterPosDegs = mean(theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs(inputConeIndices,:),1);

  % Include surround cones whose distance from the RF center are up to 
  % maxSurroundSupportFactor * surround characteristic radius of the C&K data at the cell's temporal equivalent eccentricity
  temporalEquivalentEccDegs = theMRGCMosaic.temporalEquivalentEccentricityForEccXYDegs(cachedData.RFcenterPosDegs);
  cachedData.maxSurroundSupportDegs = ...
    poolingOptimizationParamsStruct.maxSurroundSupportFactor  * ...
    RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(sqrt(sum(temporalEquivalentEccDegs.^2,2)));
  cachedData.minSurroundRadiusDegs = min(theMRGCMosaic.inputConeMosaic.coneRFspacingsDegs(:));

  % Determine the -within the max surround support- status of all cones in the cone mosaic
  coneDistancesFromRFCenter = sqrt(sum(bsxfun(@minus, theMRGCMosaic.inputConeMosaic.coneRFpositionsDegs, cachedData.RFcenterPosDegs).^2,2));
  conesAreWithinMaxSurroundSupport = (coneDistancesFromRFCenter <= cachedData.maxSurroundSupportDegs);

  % Determine the -surround connectable- status of all cones in the cone mosaic
  conesAreSurroundConnectable = false(numel(theMRGCMosaic.inputConeMosaic.coneTypes),1);
  parfor iCone = 1:numel(theMRGCMosaic.inputConeMosaic.coneTypes)
    if (ismember(theMRGCMosaic.inputConeMosaic.coneTypes(iCone), poolingOptimizationParamsStruct.coneTypesToBePooled))
      conesAreSurroundConnectable(iCone) = true;
    end
  end


  cachedData.indicesOfAllConesWithinMaxSurroundSupport = find(conesAreWithinMaxSurroundSupport);
  cachedData.indicesOfConnectableConesWithinMaxSurroundSupport = find(conesAreWithinMaxSurroundSupport & conesAreSurroundConnectable);
  cachedData.indicesOfUnconnectableConesWithinMaxSurroundSupport = find(conesAreWithinMaxSurroundSupport & ~conesAreSurroundConnectable);

  cachedData.distancesOfConnectableConesWithinMaxSurroundSupport = ...
    coneDistancesFromRFCenter(cachedData.indicesOfConnectableConesWithinMaxSurroundSupport);
  cachedData.distancesOfUnconnectableConesWithinMaxSurroundSupport = ...
    coneDistancesFromRFCenter(cachedData.indicesOfUnconnectableConesWithinMaxSurroundSupport);

end