function [theConeMosaicSTFresponses, theConeMosaicNullResponses] = ...
    computeConeMosaicSTFresponses(theConeMosaic, theOptics, ...
                                       thePresentationDisplay, ...
                                       stimParams, stimPositionDegs, ...
                                       useParfor, visualizeResponses)

    % Allocate memory
    conesNum = size(theConeMosaic.coneRFpositionsDegs,1);

    % Single precision responses
    theConeMosaicSTFresponses = zeros(...
        numel(stimParams.orientationsTested), ...
        numel(stimParams.spatialFrequenciesTested), ...
        numel(stimParams.spatialPhasesDegs), ...
        conesNum, ...
        'single');

    theConeMosaicNullResponses = [];

    if (visualizeResponses)
        useParfor = false;
    end

    if (visualizeResponses)
        hFig = figure(10);
        ax = subplot('Position', [0.1 0.1 0.85 0.85]);
    end

    % This is necessary to avoid the background being modulated when the
    % stimulus is smaller than the mosaic. This is due to the way oiCompute
    % does the padding.
    padOIwithZeros = ~true;
    if (padOIwithZeros)
        theOptics = oiSet(theOptics, 'pad', struct('value', 'zero photons'));
    end

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
                theOptics  = oiCompute(theNullStimulusScene, theOptics);

                % Compute the cone mosaic responses
                theConeMosaicNullResponses = theConeMosaic.compute(theOptics, ...
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
                theOI = theOptics;

                parfor iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)
                    % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};

                    % Compute the optical image of the frame scene
                    theCurrentOI = oiCompute(theFrameScene, theOI);

                    % Compute the cone mosaic responses
                    noiseFreeAbsorptionsCount = ...
                        theConeMosaic.compute(theCurrentOI, ...
                        'padOIwithZeros', padOIwithZeros, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

                    theFrameResponses(iFrame,:) = single(noiseFreeAbsorptionsCount(1,1,:));
                end
            else
                for iFrame = 1:numel(theCurrentStimParams.spatialPhasesDegs)

                    tic
                     % Get scene corresponding to this stimulus frame
                    theFrameScene = theDriftingGratingFrameScenes{iFrame};

                    % Compute the optical image of the frame scene
                    theOptics = oiCompute(theFrameScene, theOptics);

                    % Compute the cone mosaic responses
                    noiseFreeAbsorptionsCount = ...
                        theConeMosaic.compute(theOptics, ...
                        'padOIwithZeros', padOIwithZeros, ...
                        'opticalImagePositionDegs', stimPositionDegs, ...
                        'nTrials', 1);

                    theFrameResponses(iFrame,:) = single(noiseFreeAbsorptionsCount(1,1,:));

                    if (visualizeResponses)
                        theConeMosaicContrastResponses = ...
                            bsxfun(@times, bsxfun(@minus, noiseFreeAbsorptionsCount, theConeMosaicNullResponses), ...
                            normalizingResponses);
                        theConeMosaic.visualize(...
                            'figureHandle', hFig, ...
                            'axesHandle',ax, ...
                            'activation', theConeMosaicContrastResponses, ...
                            'activationRange', [-1 1]);
                        drawnow;
                    end

                    fprintf('Ori: %d/%d, freq: %d/%d, frame: %d of %d (%2.1f seconds)\n', ...
                        iOri, numel(stimParams.orientationsTested), iFreq, ...
                        numel(stimParams.spatialFrequenciesTested), iFrame, ...
                        numel(theCurrentStimParams.spatialPhasesDegs), toc)
                   
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
