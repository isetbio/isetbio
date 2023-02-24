function [modelConstants, retinalConePoolingParams, visualRcDegs] = computeOptimizationComponents(obj, theRGCindex)

    % Compute the visual and anatomical Rc for this RGC
    [visualRcDegs, anatomicalRcDegs, indicesOfConesPooledByTheRFcenter, weightsOfConesPooledByTheRFcenter] = ...
        computeRFcenterVisualRc(obj, theRGCindex);

    % Retinal cone pooling model constants
    modelConstants = struct();
    retinalConePoolingParams = struct();

    modelConstants.rmseWeightForRsRcResidual = obj.rmseWeightForRsRcResidual;
    modelConstants.rmseWeightForSCintSensResidual = obj.rmseWeightForSCintSensResidual;

    % The cone mosaic and the spectrally-weighted PSFs
    %modelConstants.theConeMosaic = obj.theRGCMosaic.inputConeMosaic;

    % The connectable cone types to the center and surroud
    modelConstants.surroundConnectableConeTypes = obj.retinalRFmodelParams.surroundConnectableConeTypes;
    modelConstants.centerConnectableConeTypes = obj.retinalRFmodelParams.centerConnectableConeTypes;
    modelConstants.inputConeMosaicConeTypes = obj.theRGCMosaic.inputConeMosaic.coneTypes;

    % Maximum support for the surround, in degrees, takes as twice the C&K
    % surrounds at the cell's eccentricity
    radialEccentricityForThisRGC = sqrt(sum(obj.theRGCMosaic.rgcRFpositionsDegs(theRGCindex,:).^2,2));
    
    % Surround cones up to 3 * surround characteristic radius of
    % Croner&Kaplan at this eccentricity
    modelConstants.maxSurroundSupportDegs = 3.0 * ...
                RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(radialEccentricityForThisRGC);

    % Compute the RF center position
    RFcenterPos = mean(obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);

    % Compute the distances of ALL cones in the input cone mosaic from the
    % RF center - THIS IS CACHED HERE
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);
    coneDistancesFromRFCenter = sqrt(coneDistancesFromRFCenterSquared);
    modelConstants.cachedData.surroundConeIndices = find(coneDistancesFromRFCenterSquared <= modelConstants.maxSurroundSupportDegs);
    modelConstants.cachedData.coneDistancesFromRFCenter = coneDistancesFromRFCenter(modelConstants.cachedData.surroundConeIndices);


    switch (obj.retinalRFmodelParams.retinalConePoolingModel)
        case MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex1Surround

            modelConstants.weightsComputeFunctionHandle = @MosaicPoolingOptimizer.conePoolingCoefficientsForArbitraryCenterDoubleExpSurround;
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;

            % From the 4 cells in Figure 6 of Packer & Dacey (2002)
            RnarrowToRwideRatios  = [152/515 170/718 115/902 221/1035];
            NWvolumeRatios        = [1.0     0.8     0.3     0.2];

            % Range for RwDegs, based to characteristic radius of cones and # of cones in RF center
            RwDegsInitial    = 4.00 * 1/mean(RnarrowToRwideRatios) * anatomicalRcDegs * sqrt(2.3) * sqrt(numel(modelConstants.indicesOfCenterCones));
            RwDegsLowerBound = max([0.02 0.01*RwDegsInitial]);
            RwDegsUpperBound = min([modelConstants.maxSurroundSupportDegs 4*RwDegsInitial]);

            %                                        Kc      Ks/KcRatio    narrowToWideFieldVolumeRatio  RwideDegs            RnarrowToRwideRatio
            retinalConePoolingParams.names =         {'Kc',  'KsKcRatio',  'VnVwRatio',                  'RwDegs',             'RnRwRatio'};
            retinalConePoolingParams.scaling =       {'log', 'log',        'log',                        'linear',                'log'};
            retinalConePoolingParams.initialValues = [1.       0.06        mean(NWvolumeRatios)           RwDegsInitial         mean(RnarrowToRwideRatios)];
            retinalConePoolingParams.lowerBounds   = [0.5      0.005       min(NWvolumeRatios)            RwDegsLowerBound      min(RnarrowToRwideRatios)];
            retinalConePoolingParams.upperBounds   = [2        1e0         max(NWvolumeRatios)            RwDegsUpperBound      max(RnarrowToRwideRatios)];

            H1cellIndex = 1;
            parameterTolerance = 0.3;
                
            idx = find(ismember(retinalConePoolingParams.names, 'VnVwRatio'));
            measuredValue = NWvolumeRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue - parameterTolerance;
            retinalConePoolingParams.upperBounds(idx) = measuredValue + parameterTolerance;

            idx = find(ismember(retinalConePoolingParams.names, 'RnRwRatio'));
            measuredValue = RnarrowToRwideRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue - parameterTolerance;
            retinalConePoolingParams.upperBounds(idx) = measuredValue + parameterTolerance;

        otherwise
            error('Retinal cone pooling model ''%s'' is not implemented\n', obj.retinalRFmodelParams.retinalConePoolingModel)
    end % switch (obj.retinalRFmodelParams.retinalConePoolingModel)
end

function [visualRcDegs, anatomicalRcDegs, inputConeIndices, inputConeWeights] = ...
    computeRFcenterVisualRc(obj, theTargetRGCindex)

    inputConeIndices = find(squeeze(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,theTargetRGCindex))>0.001);
    inputConeWeights = full(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(inputConeIndices, theTargetRGCindex));

    theCenterResponsesAcrossAllOrientationsAndSpatialFrequencies = sum(bsxfun(@times, ...
        obj.inputConeMosaicVisualSTFdata.responseModulations(:,:,:,inputConeIndices), ...
        reshape(inputConeWeights, [1 1 1 numel(inputConeIndices)])),4);

    % The STF of center
    theOptimalCenterSTF = obj.optimalSTFfromResponsesToAllOrientationsAndSpatialFrequencies( ...
                    theCenterResponsesAcrossAllOrientationsAndSpatialFrequencies);

    % An estimate of the anatomical RcDegs
    anatomicalRcDegs = sqrt(numel(inputConeWeights)) * ...
                       obj.theRGCMosaic.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor * ...
                       mean(obj.theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs(inputConeIndices));

    % Some initial estimate of the visual Rc (this is clearly off as it
    % does not account the effect of optics), but it is better than
    % nothing. The Gaussian STF model is only a 2-parameter model (gain,
    % Rc), so it should be able to find the Rc without getting stuct to a
    % local minimum
    initialRcDegs = anatomicalRcDegs;
    
    % Fit the visual STF with a DoG model
    fittedParamsStruct = MosaicPoolingOptimizer.fitGaussianToSubregionSTF(...
                      obj.inputConeMosaicVisualSTFdata.spatialFrequenciesTested, ...
                      theOptimalCenterSTF, ...
                      initialRcDegs, ...
                      [], ...
                      obj.multiStartsNumDoGFit);

    visualRcDegs = fittedParamsStruct.finalValues(2);
end