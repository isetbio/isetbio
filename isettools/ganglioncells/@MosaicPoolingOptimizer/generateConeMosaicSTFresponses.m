function generateConeMosaicSTFresponses(obj, gridNodeIndex, stimSizeDegs, ...
    responsesFileName, varargin)

    p = inputParser;
    p.addParameter('useParfor', false, @islogical);
    p.addParameter('visualizedResponses', false, @islogical);

    p.parse(varargin{:});
    useParfor = p.Results.useParfor;
    visualizeResponses = p.Results.visualizedResponses;

    if (isempty(gridNodeIndex))
        retinalImageResolutionDegs = retinalResolutionFromConeApertureDiameter(obj, []);
        % Stimulus centered at the RGC mosaic mosaic
        stimPositionDegs = obj.theRGCMosaic.eccentricityDegs;
        % 10% larger than the cone mosaic
        stimSizeDegs = 1.1*obj.theRGCMosaic.inputConeMosaic.sizeDegs;
        fprintf('Computing STF responses at %2.1f,%2.1f degs over a region of  %2.1f,%2.1f degs with a retinal resolution of %2.2f arc min', ...
            stimPositionDegs(1), stimPositionDegs(2), stimSizeDegs(1), stimSizeDegs(2), retinalImageResolutionDegs*60);

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
   
        fprintf('Computing STF responses at %2.1f,%2.1f degs over a region of  %2.1f,%2.1f degs with a retinal resolution of %2.2f arc min', ...
            stimPositionDegs(1), stimPositionDegs(2), stimSizeDegs(1), stimSizeDegs(2), retinalImageResolutionDegs*60);
    end


    % Generate components for running the STF mapping experiment
    [stimParams, thePresentationDisplay] = obj.setupSTFmappingExperiment(...
        stimSizeDegs, ...
        retinalImageResolutionDegs);

    % Compute cone mosaic STF responses
    [theConeMosaicSTFresponses, theConeMosaicNullResponses] = ...
        computeFocalConeMosaicSTFresponses(obj, stimParams, ...
                                           thePresentationDisplay, ...
                                           stimPositionDegs, ...
                                           useParfor, ...
                                           visualizeResponses);

    % Save computed cone mosaic STF responses to disk
    orientationsTested = stimParams.orientationsTested;
    spatialFrequenciesTested = stimParams.spatialFrequenciesTested;
    coneContrasts = stimParams.coneContrasts;
    spatialPhasesDegs = stimParams.spatialPhasesDegs;
    theNativeOpticsParams = obj.theRGCMosaic.theNativeOpticsParams;

    save(responsesFileName, 'theNativeOpticsParams', ...
        'theConeMosaicSTFresponses', 'theConeMosaicNullResponses', ...
         'orientationsTested', 'spatialFrequenciesTested', ...
         'spatialPhasesDegs', 'coneContrasts', '-v7.3');

    fprintf('Saved computed cone mosaic STF responses to %s\n', responsesFileName);

end

function [theConeMosaicSTFresponses, theConeMosaicNullResponses] = ...
    computeFocalConeMosaicSTFresponses(obj, stimParams, ...
                                       thePresentationDisplay, stimPositionDegs, ...
                                       useParfor, visualizeResponses)

    % Allocate memory
    conesNum = size(obj.theRGCMosaic.inputConeMosaic.coneRFpositionsDegs,1);

    % Single precision responses
    theConeMosaicSTFresponses = zeros(...
        numel(stimParams.orientationsTested), ...
        numel(stimParams.spatialFrequenciesTested), ...
        numel(stimParams.spatialPhasesDegs), ...
        conesNum, ...
        'single');

    theConeMosaicNullResponses = [];

    if (useParfor)
        visualizeResponses = false;
    end

    if (visualizeResponses)
        hFig = figure(10);
        ax = subplot('Position', [0.1 0.1 0.85 0.85]);
    end

    % Retrieve the native optics
    nativeOI = obj.theRGCMosaic.theNativeOptics;

    % This is necessary to avoid the background being modulated when the
    % stimulus is smaller than the mosaic. This is due to the way oiCompute
    % does the padding.
    padOIwithZeros = true;
    if (padOIwithZeros)
        nativeOI = oiSet(nativeOI, 'pad', struct('value', 'zero photons'));
    end

    % Retrieve the input cone mosaic
    inputConeMosaic = obj.theRGCMosaic.inputConeMosaic;

    % Go through all stimulus orientations
    for iOri = 1:numel(stimParams.orientationsTested)
        stimParams.orientationDegs = stimParams.orientationsTested(iOri);

        fprintf('Computing cone mosaic STFs for the %d degs orientation patterns.\n', ...
                stimParams.orientationDegs);

        for iFreq = numel(stimParams.spatialFrequenciesTested):-1:1

            theCurrentStimParams = stimParams;
            theCurrentStimParams.spatialFrequencyCPD = stimParams.spatialFrequenciesTested(iFreq);

            % Generate spatial modulation patterns for each stimulus frame
            theDriftingGratingSpatialModulationPatterns = rfMappingStimulusGenerator.driftingGratingFrames(theCurrentStimParams);

            % Generate scenes for the different spatial phases
            [theDriftingGratingFrameScenes, theNullStimulusScene] = ...
                rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                    thePresentationDisplay, theCurrentStimParams, theDriftingGratingSpatialModulationPatterns, ...
                    'validateScenes', false);

            % Compute the cone mosaic responses to the null scene OI
            if (isempty(theConeMosaicNullResponses))
                % Compute the optical image of the null scene
                nativeOI  = oiCompute(theNullStimulusScene, nativeOI);

                % Compute the cone mosaic responses
                theConeMosaicNullResponses = inputConeMosaic.compute(nativeOI, ...
                    'padOIwithZeros', padOIwithZeros, ...
                    'opticalImagePositionDegs', stimPositionDegs, ...
                    'nTrials', 1);

                if (visualizeResponses)
                    coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponses== 0);
                    normalizingResponses = 1./theConeMosaicNullResponses;
                    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
                    normalizingResponses = reshape(normalizingResponses, [1 1 numel(normalizingResponses)]);
                end

            end

            % Allocate memory
            theFrameResponses = zeros(numel(theCurrentStimParams.spatialPhasesDegs), conesNum, 'single');

            % Compute the input cone mosaic responses
            if (useParfor)
                theOI = nativeOI;
                parfor iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};

                    % Compute the optical image of the frame scene
                    theCurrentOI = oiCompute(theFrameScene, theOI);

                    % Compute the cone mosaic responses
                    noiseFreeAbsorptionsCount = ...
                        inputConeMosaic.compute(theCurrentOI, ...
                        'padOIwithZeros', padOIwithZeros, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

                    theFrameResponses(iFrame,:) = single(noiseFreeAbsorptionsCount(1,1,:));
                end
            else
                for iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};

                    % Compute the optical image of the frame scene
                    nativeOI = oiCompute(theFrameScene, nativeOI);

                    % Compute the cone mosaic responses
                    noiseFreeAbsorptionsCount = ...
                        inputConeMosaic.compute(nativeOI, ...
                        'padOIwithZeros', padOIwithZeros, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

                    theFrameResponses(iFrame,:) = single(noiseFreeAbsorptionsCount(1,1,:));

                    if (visualizeResponses)
                        theConeMosaicContrastResponses = ...
                            bsxfun(@times, bsxfun(@minus, noiseFreeAbsorptionsCount, theConeMosaicNullResponses), ...
                            normalizingResponses);
                        inputConeMosaic.visualize(...
                            'figureHandle', hFig, ...
                            'axesHandle',ax, ...
                            'activation', theConeMosaicContrastResponses, ...
                            'activationRange', [-1 1]);
                        drawnow;
                    end
                end  % iFrame
            end


            % Save theFrameResponses
            theConeMosaicSTFresponses(iOri, iFreq,:,:) = theFrameResponses;

            % Save memory
            theDriftingGratingFrameScenes = [];
            theNullStimulusScene = [];

        end % iFreq
    end % iOri
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
