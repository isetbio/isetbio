function generateInputConeMosaicSubspaceRFmappingLinearResponses(theRGCMosaic, ...
    opticsToEmploy, maxSFcyclesPerDegree, stimSizeDegs, stimXYpositionDegs, responsesFileName, varargin)

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

    % Retrieve the optics to employ
    switch (opticsToEmploy)
        case {'native', 'adaptive optics'}
             % Retrieve the native optics
             theOptics = theRGCMosaic.theNativeOptics;
             
        case 'custom'
            % Retrieve the custom optics
             theOptics = theRGCMosaic.theCustomOptics;

        otherwise
             error('Unknown optics: ''%s''.', opticsToEmploy);
    end

    [theConeMosaicSubspaceLinearResponses, theConeMosaicNullResponses, ...
          HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices] = ...
          MosaicPoolingOptimizer.computeConeMosaicSubspaceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimXYpositionDegs, ...
                                           'visualizeResponses', visualizedResponses, ...
                                           'parPoolSize', parPoolSize);

    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);

    switch (opticsToEmploy)
        case {'native', 'adaptive optics'}
             theNativeOpticsParams = theRGCMosaic.theNativeOpticsParams;
             save(responsesFileName, 'theNativeOpticsParams', ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearResponses', 'theConeMosaicNullResponses', ...
                 '-v7.3');

        case 'custom'
             theCustomOpticsParams = theRGCMosaic.theCustomOpticsParams;
             save(responsesFileName, 'theCustomOpticsParams', ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearResponses', 'theConeMosaicNullResponses', ...
                 '-v7.3');
        otherwise
             error('Unknown optics: ''%s''.', opticsToEmploy);
    end

    fprintf('Saved computed cone mosaic SUBSPACE RF mapping linear responses to %s\n', responsesFileName);

end
