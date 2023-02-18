function theSTFdata = visualSTFfromCronerKaplanAnalysisOfVisualRF(obj, theVisualRF)

    % Rotate theVisualRFmap according to the rotation
    % that maximizes horizontal resolution of the targetVisualRFmap
    theRotatedVisualRF = RTVF.bestHorizontalResolutionRFmap(...
        theVisualRF, obj.bestHorizontalResolutionRotationDegs);

    % Integrate along Y to generate the X-axis line weighting function
    spatialSupportDegsY = obj.spectrallyWeightedPSFData.spatialSupportForRFmapYdegs;
    spatialSampleSize = spatialSupportDegsY(2)-spatialSupportDegsY(1);
    theLineWeightingFunctionX = sum(theRotatedVisualRF,1) * spatialSampleSize;

    % Compute the visual STF corresponding to the X-axis line weighting function
    spatialSupportDegsX = obj.spectrallyWeightedPSFData.spatialSupportForRFmapXdegs;
    [theSpatialFrequencySupport, theVisualSTF] = ...
        RTVF.spatialTransferFunction(spatialSupportDegsX, theLineWeightingFunctionX);


    % Fit the visual STF with a DoG model
    [theFittedDoGModelParams, theFittedDoGModelToTheVisualSTF] = RTVF.fitDoGmodelToMeasuredSTF(...
                      theSpatialFrequencySupport, ...
                      theVisualSTF, ...
                      obj.visualRFcenterRcDegs, ...
                      obj.visualRFcenterRcDegs*[0.5 1 1.5], ...
                      obj.multiStartsNumDoGFit);


     theSTFdata = struct;
     theSTFdata.spatialFrequencySupport = theSpatialFrequencySupport;
     theSTFdata.visualSTF = theVisualSTF;
     theSTFdata.fittedDoGModelParams = theFittedDoGModelParams;
     theSTFdata.fittedDoGModelToVisualSTF = theFittedDoGModelToTheVisualSTF;
     theSTFdata.fittedRcDegs = theFittedDoGModelParams.finalValues(4);
     theSTFdata.fittedDoGModelRsRcRatio = theFittedDoGModelParams.finalValues(3);
     theSTFdata.fittedDoGModelSCIntSensRatio = theFittedDoGModelParams.finalValues(2) * (theFittedDoGModelParams.finalValues(3))^2;

end
