function generateInputConeMosaicSubspaceRFmappingResponses(theRGCMosaic, responsesFileName)

    % Determine optimal stimulus resolution so that cone aperture blur will
    % have an observable effect
    retinalImageResolutionDegs = min(theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs)/9;
    
    % Stimulus centered at the RGC mosaic mosaic
    stimPositionDegs = theRGCMosaic.eccentricityDegs;

    % 10% larger than the cone mosaic
    stimSizeDegs = 1.1*theRGCMosaic.inputConeMosaic.sizeDegs;
    fprintf('Computing STF responses for stimuli positioned at %2.1f,%2.1f degs\nover a region of  %2.1f,%2.1f degs\nwith a retinal resolution of %2.2f arc min\n', ...
            stimPositionDegs(1), stimPositionDegs(2), stimSizeDegs(1), stimSizeDegs(2), retinalImageResolutionDegs*60);


    % Generate components for running the Subspace mapping experiment
    [stimParams, thePresentationDisplay] = MosaicPoolingOptimizer.setupSubspaceRFmappingExperiment(...
        theRGCMosaic.inputConeMosaic, ...
        stimSizeDegs, retinalImageResolutionDegs);

    % Retrieve the native optics
    theOptics = theRGCMosaic.theNativeOptics;

    % Retrieve the input cone mosaic
    theConeMosaic = theRGCMosaic.inputConeMosaic;

    % Compute the subspace responses of the input cone mosaic
    useParfor = true;
    visualizeResponses = false;
    [theConeMosaicSubspaceResponses, theConeMosaicNullResponses, ...
        HartleySpatialModulationPatterns, lIndices, mIndices] = ...
        MosaicPoolingOptimizer.computeConeMosaicSubspaceRFmappingResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimPositionDegs, ...
                                           useParfor, ...
                                           visualizeResponses);

    theNativeOpticsParams = theRGCMosaic.theNativeOpticsParams;

    save(responsesFileName, 'theNativeOpticsParams', ...
                'HartleySpatialModulationPatterns', 'lIndices', 'mIndices', ...
                'theConeMosaicSubspaceResponses', 'theConeMosaicNullResponses', ...
                 '-v7.3');

    fprintf('Saved computed cone mosaic SUBSPACE RF mapping responses to %s\n', responsesFileName);

end
