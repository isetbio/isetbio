function generateInputConeMosaicSubspaceRFmappingLinearResponses(theRGCMosaic, ...
    theOptics, maxSFcyclesPerDegree, stimSizeDegs, stimXYpositionDegs, responsesFileName, varargin)

    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizedResponses = p.Results.visualizedResponses;

    % Generate components for running the Subspace mapping experiment
    wavelengthSupport = theRGCMosaic.inputConeMosaic.wave;
   
    [stimParams, thePresentationDisplay] = MosaicPoolingOptimizer.setupSubspaceRFmappingExperiment(...
        wavelengthSupport, stimSizeDegs, maxSFcyclesPerDegree);

    % Retrieve the input cone mosaic
    theConeMosaic = theRGCMosaic.inputConeMosaic;

    % Compute the subspace responses of the input cone mosaic.
    % Note: These are linear responses, i.e., 
    % Response(forward polarity Hartley) - Response (inverse polarity Hartley)

    [theConeMosaicSubspaceLinearResponses, theConeMosaicSubspaceEnergyResponses, theConeMosaicNullResponses, ...
          HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices] = ...
          MosaicPoolingOptimizer.computeConeMosaicSubspaceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimXYpositionDegs, ...
                                           'visualizeResponses', visualizedResponses, ...
                                           'parPoolSize', parPoolSize);


    theNativeOpticsParams = theRGCMosaic.theNativeOpticsParams;
    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);
    
    save(responsesFileName, 'theNativeOpticsParams', ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearResponses', 'theConeMosaicSubspaceEnergyResponses', 'theConeMosaicNullResponses', ...
                 '-v7.3');

    fprintf('Saved computed cone mosaic SUBSPACE RF mapping linear responses to %s\n', responsesFileName);

end
