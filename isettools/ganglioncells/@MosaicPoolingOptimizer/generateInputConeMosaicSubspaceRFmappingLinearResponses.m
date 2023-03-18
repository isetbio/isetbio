function generateInputConeMosaicSubspaceRFmappingLinearResponses(theRGCMosaic, responsesFileName)

    % Generate components for running the Subspace mapping experiment
    maxSFcyclesPerDegree = 30;
    wavelengthSupport = theRGCMosaic.inputConeMosaic.wave;
    [stimParams, thePresentationDisplay] = MosaicPoolingOptimizer.setupSubspaceRFmappingExperiment(...
        wavelengthSupport, stimSizeDegs, maxSFcyclesPerDegree);

    % Retrieve the native optics
    theOptics = theRGCMosaic.theNativeOptics;

    % Retrieve the input cone mosaic
    theConeMosaic = theRGCMosaic.inputConeMosaic;

    % Compute the subspace responses of the input cone mosaic.
    % Note: These are linear responses, i.e., 
    % Response(forward polarity Hartley) - Response (inverse polarity Hartley)

    [theConeMosaicSubspaceResponses, theConeMosaicNullResponses, ...
          HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices] = ...
          MosaicPoolingOptimizer.computeConeMosaicSubspaceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           theRGCMosaic.eccentricityDegs, ...
                                           useParfor, ...
                                           'parPoolSize', []);


    theNativeOpticsParams = theRGCMosaic.theNativeOpticsParams;
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    
    save(responsesFileName, 'theNativeOpticsParams', ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceResponses', 'theConeMosaicNullResponses', ...
                 '-v7.3');

    fprintf('Saved computed cone mosaic SUBSPACE RF mapping linear responses to %s\n', responsesFileName);

end
