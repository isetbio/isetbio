function generateInputConeMosaicSTFresponses(obj, gridNodeIndex, stimSizeDegs, ...
    responsesFileName, varargin)

    p = inputParser;
    p.addParameter('useParfor', false, @islogical);
    p.addParameter('visualizedResponses', false, @islogical);
    p.addParameter('opticsToEmploy', 'native', @(x)(ismember(x, {'native', 'custom'})));

    p.parse(varargin{:});
    useParfor = p.Results.useParfor;
    visualizeResponses = p.Results.visualizedResponses;
    opticsToEmploy = p.Results.opticsToEmploy;

    if (isempty(gridNodeIndex))
        % Determine optimal stimulus resolution so that cone aperture blur will
        % have an observable effect
        retinalImageResolutionDegs = retinalResolutionFromConeApertureDiameter(obj, []);
        % Stimulus centered at the RGC mosaic mosaic
        stimPositionDegs = obj.theRGCMosaic.eccentricityDegs;
        % 10% larger than the cone mosaic
        stimSizeDegs = 1.1*obj.theRGCMosaic.inputConeMosaic.sizeDegs;
        fprintf('Computing STF responses for stimuli positioned at %2.1f,%2.1f degs\nover a region of  %2.1f,%2.1f degs\nwith a retinal resolution of %2.2f arc min\nand %s optics', ...
            stimPositionDegs(1), stimPositionDegs(2), stimSizeDegs(1), stimSizeDegs(2), retinalImageResolutionDegs*60, opticsToEmploy);

    else
        % Retrieve the RGC indices for the L-center and M-center RGC at the desired grid node
        targetRGCindices = [...
            obj.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex) ...
            obj.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex)];
    
        % Determine position of node (mean of L-center and M-center RGC)
        stimPositionDegs = mean(obj.theRGCMosaic.rgcRFpositionsDegs(targetRGCindices,:),1);
        
        % Determine optimal stimulus resolution so that cone aperture blur will
        % have an observable effect
        retinalImageResolutionDegs = retinalResolutionFromConeApertureDiameter(obj, targetRGCindices);
   
        fprintf('Computing STF responses for stimuli positioned at %2.1f,%2.1f degs\nover a region of  %2.1f,%2.1f degs\nwith a retinal resolution of %2.2f arc min\n and % optics', ...
            stimPositionDegs(1), stimPositionDegs(2), stimSizeDegs(1), stimSizeDegs(2), retinalImageResolutionDegs*60, opticsToEmploy);
    end


    % Generate components for running the STF mapping experiment
    [stimParams, thePresentationDisplay] = MosaicPoolingOptimizer.setupSTFmappingExperiment(...
        obj.theRGCMosaic.inputConeMosaic, ...
        stimSizeDegs, ...
        retinalImageResolutionDegs);

    % Compute cone mosaic STF responses
    switch (opticsToEmploy)
        case 'native'
             % Retrieve the native optics
             theOptics = obj.theRGCMosaic.theNativeOptics;
             
        case 'custom'
            % Retrieve the custom optics
             theOptics = obj.theRGCMosaic.theCustomOptics;

        otherwise
             error('Unknown optics: ''%s''.', opticsToEmploy);
    end
    fprintf('Will save responses in %s.\n', responsesFileName);

   
    % Retrieve the input cone mosaic
    theConeMosaic = obj.theRGCMosaic.inputConeMosaic;

    [theConeMosaicSTFresponses, theConeMosaicNullResponses] = ...
        MosaicPoolingOptimizer.computeConeMosaicSTFresponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimPositionDegs, ...
                                           useParfor, ...
                                           visualizeResponses);

    % Save computed cone mosaic STF responses to disk
    orientationsTested = stimParams.orientationsTested;
    spatialFrequenciesTested = stimParams.spatialFrequenciesTested;
    coneContrasts = stimParams.coneContrasts;
    spatialPhasesDegs = stimParams.spatialPhasesDegs;


    switch (opticsToEmploy)
        case 'native'
             theNativeOpticsParams = obj.theRGCMosaic.theNativeOpticsParams;
             save(responsesFileName, 'theNativeOpticsParams', ...
                'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
                 'orientationsTested', 'spatialFrequenciesTested', ...
                 'spatialPhasesDegs', 'coneContrasts', '-v7.3');
        case 'custom'
             theCustomOpticsParams = obj.theRGCMosaic.theCustomOpticsParams;
             save(responsesFileName, 'theCustomOpticsParams', ...
                'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
                 'orientationsTested', 'spatialFrequenciesTested', ...
                 'spatialPhasesDegs', 'coneContrasts', '-v7.3');
        otherwise
             error('Unknown optics: ''%s''.', opticsToEmploy);
    end

    fprintf('Saved computed cone mosaic STF responses using **%s** optics to %s\n', opticsToEmploy, responsesFileName);

end


function retinalImageResolutionDegs = retinalResolutionFromConeApertureDiameter(obj, targetRGCindices)

    if (isempty(targetRGCindices))
        retinalImageResolutionDegs = min(obj.theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs)/9;
    else
        % Find cone aperture size of the input cones
        coneIndices = find(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,targetRGCindices(1)) > 0.001);
        coneIndices = cat(1, coneIndices, ...
                      find(obj.theRGCMosaic.rgcRFcenterConeConnectivityMatrix(:,targetRGCindices(2)) > 0.001));
        
        retinalImageResolutionDegs = mean(obj.theRGCMosaic.inputConeMosaic.coneApertureDiametersDegs(coneIndices))/13; 
    end

end