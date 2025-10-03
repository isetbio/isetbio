function generateInputConeMosaicSubspaceRFmappingLinearResponses(theRGCMosaic, ...
    opticsToEmploy, stimSizeDegs, stimXYpositionDegs, ...
    stimulusChromaticity, coneFundamentalsOptimizedForStimPosition, ...
    responsesFileName, varargin)

    p = inputParser;
    p.addParameter('maxSFLimit', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('rfMappingPixelMagnificationFactor', 1, @(x)(isscalar(x)&&(x>=1)));
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizedResponses', false, @islogical);
    p.parse(varargin{:});

    maxSFLimit = p.Results.maxSFLimit;
    rfMappingPixelMagnificationFactor = p.Results.rfMappingPixelMagnificationFactor;
    parPoolSize = p.Results.parPoolSize;
    visualizedResponses = p.Results.visualizedResponses;

    % Generate components for running the Subspace mapping experiment
    wavelengthSupport = theRGCMosaic.inputConeMosaic.wave;
   

    % Find RGCs within the stimulus region
    if (numel(stimSizeDegs) == 2)
        widthDegs = stimSizeDegs(1);
        heightDegs = stimSizeDegs(2);
    else
        widthDegs = stimSizeDegs(1);
        heightDegs = stimSizeDegs(1);
    end

    theStimulusRegion = regionOfInterest(...
        'geometryStruct', struct(...
            'units', 'degs', ...
            'shape', 'rect', ...
            'center', stimXYpositionDegs, ...
            'width', widthDegs, ...
            'height', heightDegs , ...
            'rotation', 0.0...
        ));

    targetRGCindices = theStimulusRegion.indicesOfPointsInside(theRGCMosaic.rgcRFpositionsDegs);

    % Determine optimal stimulus resolution so that cone aperture blur will
    % have an observable effect
    optimalRetinalPixelSizeDegs = MosaicPoolingOptimizer.retinalResolutionFromConeApertureDiameter(theRGCMosaic, targetRGCindices);

    employedRetinalPixelSizeDegs = rfMappingPixelMagnificationFactor * optimalRetinalPixelSizeDegs;

    [stimParams, thePresentationDisplay] = MosaicPoolingOptimizer.setupSubspaceRFmappingExperiment(...
        wavelengthSupport, stimSizeDegs, employedRetinalPixelSizeDegs, maxSFLimit, stimulusChromaticity);

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

    if (coneFundamentalsOptimizedForStimPosition)
        % Compute custom cone fundamentals
        maxConesNumForAveraging = 3;
        customConeFundamentals = MosaicPoolingOptimizer.coneFundamentalsAtTargetPositionWithinConeMosaic(...
            theConeMosaic, theOptics, stimXYpositionDegs, stimParams.stimSizeDegs, maxConesNumForAveraging);
    else
        customConeFundamentals = [];
    end

    [theConeMosaicSubspaceLinearModulationResponses, theConeMosaicNullResponses, ...
          HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices] = ...
          MosaicPoolingOptimizer.computeConeMosaicSubspaceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimXYpositionDegs, ...
                                           'customConeFundamentals', customConeFundamentals, ...
                                           'visualizeResponses', visualizedResponses, ...
                                           'parPoolSize', parPoolSize);

    HartleySpatialModulationPatterns = single(HartleySpatialModulationPatterns);

    switch (opticsToEmploy)
        case {'native', 'adaptive optics'}
             theNativeOpticsParams = theRGCMosaic.theNativeOpticsParams;
             save(responsesFileName, 'theNativeOpticsParams', ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearModulationResponses', 'theConeMosaicNullResponses', ...
                 '-v7.3');

        case 'custom'
             theCustomOpticsParams = theRGCMosaic.theCustomOpticsParams;
             save(responsesFileName, 'theCustomOpticsParams', ...
                'HartleySpatialModulationPatterns', 'spatialSupportDegs', 'stimParams', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceLinearModulationResponses', 'theConeMosaicNullResponses', ...
                 '-v7.3');
        otherwise
             error('Unknown optics: ''%s''.', opticsToEmploy);
    end

    fprintf('Saved computed cone mosaic SUBSPACE RF mapping linear responses to %s\n', responsesFileName);

end
