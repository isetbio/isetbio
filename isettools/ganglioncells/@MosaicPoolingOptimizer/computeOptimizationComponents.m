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
    
    % Surround cones up to 2 * surround characteristic radius of
    % Croner&Kaplan at this eccentricity
    modelConstants.maxSurroundSupportDegs = MosaicPoolingOptimizer.maxSurroundSupportFactor  * ...
                RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(radialEccentricityForThisRGC);

    % Compute the RF center position
    RFcenterPos = mean(obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);

    % Compute the distances of ALL cones in the input cone mosaic from the
    % RF center - THIS IS CACHED HERE
    coneDistancesFromRFCenterSquared = sum(bsxfun(@minus, obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs, RFcenterPos).^2,2);
    coneDistancesFromRFCenter = sqrt(coneDistancesFromRFCenterSquared);
    modelConstants.cachedData.surroundConeIndices = find(coneDistancesFromRFCenter <= modelConstants.maxSurroundSupportDegs);
    modelConstants.cachedData.coneDistancesFromRFCenter = coneDistancesFromRFCenter(modelConstants.cachedData.surroundConeIndices);


    switch (obj.retinalRFmodelParams.conePoolingModel)
        case {
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex1Surround, ...
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex2Surround, ...
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex3Surround, ...
                MosaicPoolingOptimizer.ArbitraryCenter_DoubleExpH1cellIndex4Surround, ...
                }

            % From the 4 cells in Figure 6 of Packer & Dacey (2002)
            RnarrowToRwideRatios  = MosaicPoolingOptimizer.PackerDacey2002_H1params.RnarrowToRwideRatios;
            NWvolumeRatios        = MosaicPoolingOptimizer.PackerDacey2002_H1params.NWvolumeRatios;

            modelConstants.retinalConePoolingModel = obj.retinalRFmodelParams.conePoolingModel;
            modelConstants.weightsComputeFunctionHandle = @MosaicPoolingOptimizer.conePoolingCoefficientsForArbitraryCenterDoubleExpSurround;
            modelConstants.indicesOfCenterCones = indicesOfConesPooledByTheRFcenter;
            modelConstants.weightsOfCenterCones = weightsOfConesPooledByTheRFcenter;

            
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

            if (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex1'))
                H1cellIndex = 1;
                fprintf('Fitting using params from the FIRST H1 cell of Dacey and Packer\n');
            elseif (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex2'))
                H1cellIndex = 2;
                fprintf('Fitting using params from the SECOND H1 cell of Dacey and Packer\n');
            elseif (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex3'))
                H1cellIndex = 3;
                fprintf('Fitting using params from the THIRD H1 cell of Dacey and Packer\n');
            elseif (contains(modelConstants.retinalConePoolingModel, 'H1cellIndex4'))
                H1cellIndex = 4;
                fprintf('Fitting using params from the FOURTH H1 cell of Dacey and Packer\n');
            else
                error('Could not determine H1cellIndex from the retinal cone pooling model: ''%s''.\n', ...
                    modelConstants.retinalConePoolingModel);
            end

            idx = find(ismember(retinalConePoolingParams.names, 'VnVwRatio'));
            measuredValue = NWvolumeRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue * (1 - obj.retinalRFmodelParams.H1parameterTolerance);
            retinalConePoolingParams.upperBounds(idx) = measuredValue * (1 + obj.retinalRFmodelParams.H1parameterTolerance);

            idx = find(ismember(retinalConePoolingParams.names, 'RnRwRatio'));
            measuredValue = RnarrowToRwideRatios(H1cellIndex);
            retinalConePoolingParams.initialValues(idx) = measuredValue;
            retinalConePoolingParams.lowerBounds(idx) = measuredValue * (1 - obj.retinalRFmodelParams.H1parameterTolerance);
            retinalConePoolingParams.upperBounds(idx) = measuredValue * (1 + obj.retinalRFmodelParams.H1parameterTolerance);
           

        otherwise
            error('Retinal cone pooling model ''%s'' is not implemented\n', obj.retinalRFmodelParams.conePoolingModel)
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