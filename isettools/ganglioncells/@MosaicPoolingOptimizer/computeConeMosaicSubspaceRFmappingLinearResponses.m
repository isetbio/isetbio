function [theConeMosaicSubspaceLinearResponses, theConeMosaicNullResponses, ...
    HartleySpatialModulationPatterns, spatialSupportDegs, lIndices, mIndices] = ...
    computeConeMosaicSubspaceRFmappingLinearResponses(theConeMosaic, theOptics,  ...
                                           thePresentationDisplay, ...
                                           stimParams, ...
                                           stimPositionDegs, ...
                                           varargin)

    p = inputParser;
    p.addParameter('parPoolSize', [], @(x)(isempty(x)||(isscalar(x))));
    p.addParameter('visualizeResponses', false, @islogical);
    p.parse(varargin{:});
    parPoolSize = p.Results.parPoolSize;
    visualizeResponses = p.Results.visualizeResponses;

    % Compute the Hartley spatial patterns
    visualizePatterns = ~true;

    if (visualizePatterns)
        % Compute spatial modulation patterns for the Hartley set
        [HartleySpatialModulationPatterns, lIndices, mIndices] = ...
                rfMappingStimulusGenerator.HartleyModulationPatterns(...
                    stimParams.omega, stimParams.stimSizeDegs, stimParams.pixelSizeDegs, ...
                    'parPoolSize', 0, ...
                    'visualizePatterns', true);
    else
        [HartleySpatialModulationPatterns, lIndices, mIndices] = ...
            rfMappingStimulusGenerator.HartleyModulationPatterns(...
                stimParams.omega, stimParams.stimSizeDegs, stimParams.pixelSizeDegs);
    end

    % Compute memory requirement
    a = single(1);
    s = whos('a');
    nStim = size(HartleySpatialModulationPatterns,1);
    memRequirementGBytes = (nStim * theConeMosaic.conesNum * s.bytes)/1024/1024/1024;
    fprintf('The cone mosaic responses to the %d Hartley patterns require %2.1f GBytes of memory.\n', ...
         nStim, memRequirementGBytes)


    % Compute the null stimulus. Make sure it covers all of the cone mosaic
    fprintf('Computing null scene\n');
    nullStimParams = stimParams;
    nullStimParams.sizeDegs = 1.1*max(theConeMosaic.sizeDegs);
    [~, theNullStimulusScene, spatialSupportDegs] = rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                thePresentationDisplay, stimParams, HartleySpatialModulationPatterns, ...
                'validateScenes', false, ...
                'sceneIndexToCompute', 0);

    fprintf('Computing null scene response\n');

    % Compute the optical image of the null scene
    theOptics  = oiCompute(theNullStimulusScene, theOptics);

    % Compute the cone mosaic response to the null stimulus
    theConeMosaicNullResponses = theConeMosaic.compute(theOptics, ...
                    'opticalImagePositionDegs', theConeMosaic.eccentricityDegs, ...
                    'nTrials', 1);

    coneIndicesWithZeroNullResponse = find(theConeMosaicNullResponses== 0);
    normalizingResponses = 1./theConeMosaicNullResponses;
    normalizingResponses(coneIndicesWithZeroNullResponse) = 0;
    normalizingResponses = reshape(normalizingResponses, [1 1 numel(normalizingResponses)]);
             

    % Compute the input cone mosaic responses
    theConeMosaicSubspaceLinearResponses = zeros(nStim, theConeMosaic.conesNum, 'single');

    if ((~isempty(parPoolSize)) && (parPoolSize>1)) || (isempty(parPoolSize)) || (visualizeResponses)
         % Reset parpool
         [shutdownParPoolOnceCompleted, numWorkers] = MosaicPoolingOptimizer.resetParPool(parPoolSize);

         theOI = theOptics;
         parfor iFrame = 1:nStim
             % Generate scenes for the Hartley patterns
            fprintf('Computing cone mosaic response for Hartley pattern %d of %d (using %d parallel processes).\n', ...
                iFrame, nStim, numWorkers);

             theForwardPolarityRFMappingStimulusScenes =  rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                thePresentationDisplay, stimParams, HartleySpatialModulationPatterns, ...
                'validateScenes', false, ...
                'sceneIndexToCompute', iFrame);

             theInversePolarityRFMappingStimulusScenes = rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                thePresentationDisplay, stimParams, -HartleySpatialModulationPatterns, ...
                'validateScenes', false, ...
                'sceneIndexToCompute', iFrame);

             % Get scene corresponding to the forward polarity of the stimulus frame
             theFrameScene = theForwardPolarityRFMappingStimulusScenes{1};

             % Compute the optical image of the frame scene
             theCurrentOI = oiCompute(theFrameScene, theOI);

             
             % Compute the cone mosaic responses
             noiseFreeAbsorptionsCountForwardPolarity = ...
                        theConeMosaic.compute(theCurrentOI, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

             % Convert to modulation responses
             noiseFreeAbsorptionsCountForwardPolarity = ...
                (noiseFreeAbsorptionsCountForwardPolarity - theConeMosaicNullResponses) .* ...
                normalizingResponses;

             % Get scene corresponding to the inverse polarity of this stimulus frame
             theFrameScene = theInversePolarityRFMappingStimulusScenes{1};

             % Compute the optical image of the frame scene
             theCurrentOI = oiCompute(theFrameScene, theOI);

             % Compute the cone mosaic responses
             noiseFreeAbsorptionsCountInversePolarity = ...
                        theConeMosaic.compute(theCurrentOI, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

             % Convert to modulation responses
             noiseFreeAbsorptionsCountInversePolarity = ...
                (noiseFreeAbsorptionsCountInversePolarity - theConeMosaicNullResponses) .* ...
                normalizingResponses;


             % The linear response: forward - inverse
             theConeMosaicSubspaceLinearResponses(iFrame,:) = single(...
                 noiseFreeAbsorptionsCountForwardPolarity(1,1,:) - ...
                 noiseFreeAbsorptionsCountInversePolarity(1,1,:));


         end % iFrame

         if (shutdownParPoolOnceCompleted)
            poolobj = gcp('nocreate'); 
            delete(poolobj);
         end
         
    else
        for iFrame = 1:nStim

            fprintf('Computing cone mosaic response for Hartley pattern %d of %d (serially).\n', iFrame, nStim);

            theForwardPolarityRFMappingStimulusScenes =  rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                thePresentationDisplay, stimParams, HartleySpatialModulationPatterns, ...
                'validateScenes', false, ...
                'sceneIndexToCompute', iFrame);

            theInversePolarityRFMappingStimulusScenes = rfMappingStimulusGenerator.generateStimulusFramesOnPresentationDisplay(...
                thePresentationDisplay, stimParams, -HartleySpatialModulationPatterns, ...
                'validateScenes', false, ...
                'sceneIndexToCompute', iFrame);

            % Get scene corresponding to this forward polarity of the stimulus frame
            theFrameScene = theForwardPolarityRFMappingStimulusScenes{1};


            % Compute the optical image of the frame scene
            theOptics = oiCompute(theFrameScene, theOptics);

            % Compute the cone mosaic responses
            noiseFreeAbsorptionsCountForwardPolarity = ...
                        theConeMosaic.compute(theOptics, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

            % Convert to modulation responses
            noiseFreeAbsorptionsCountForwardPolarity = ...
                (noiseFreeAbsorptionsCountForwardPolarity - theConeMosaicNullResponses) ./ ...
                normalizingResponses;

            % Get scene corresponding to the inverse polarity of this stimulus frame
            theFrameScene = theInversePolarityRFMappingStimulusScenes{1};

            % Compute the optical image of the frame scene
            theOptics = oiCompute(theFrameScene, theOptics);

            % Compute the cone mosaic responses
            noiseFreeAbsorptionsCountInversePolarity = ...
                        theConeMosaic.compute(theOptics, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

            % Convert to modulation responses
            noiseFreeAbsorptionsCountInversePolarity = ...
                (noiseFreeAbsorptionsCountInversePolarity - theConeMosaicNullResponses) ./ ...
                normalizingResponses;

            if (visualizeResponses)
                hFig = figure(2); clf;
                set(hFig, 'Position', [10 10 2000 1200]);
                ax1 = subplot(2,1,1);
                ax2 = subplot(2,1,2);
    
                theConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', ax1, ...
                    'activation', noiseFreeAbsorptionsCountForwardPolarity);
                theConeMosaic.visualize(...
                    'figureHandle', hFig, ...
                    'axesHandle', ax2, ...
                    'activation', noiseFreeAbsorptionsCountInversePolarity);
                drawnow;
            end

            % The linear response: forward - inverse
            theConeMosaicSubspaceLinearResponses(iFrame,:) = single(...
                 noiseFreeAbsorptionsCountForwardPolarity(1,1,:) - ...
                 noiseFreeAbsorptionsCountInversePolarity(1,1,:));
        end
    end
end
