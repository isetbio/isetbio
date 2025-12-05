%
% RGCMosaicConstructor.helper.simulateExperiment.temporalTransferFunction
%
%
function temporalTransferFunction(theMRGCMosaic, theOI, ...
    TTFparamsStruct, theInputConeMosaicTTFResponsesFullFileName, ...
    theMRGCMosaicTTFResponsesFullFileName, varargin)

    p = inputParser;
    p.addParameter('computeInputConeMosaicResponses', false, @islogical);
    p.addParameter('computeInputConeMosaicResponsesBasedOnConeExcitations', true, @islogical);
    p.addParameter('computeInputConeMosaicResponsesBasedOnPhotocurrents', true, @islogical);
    p.addParameter('inspectInputConeMosaicResponses', false, @islogical);
    p.addParameter('computeMRGCMosaicResponses', false, @islogical);
    p.addParameter('debugInputConeMosaicPcurrentResponse', false, @islogical);
    p.addParameter('visualizeResponse', false, @islogical);
    p.addParameter('visualizeStimulusSequence', false, @islogical);
    p.addParameter('visualizeStimulusOnMosaic', false, @islogical);
    p.addParameter('mRGCNonLinearityParams', [], @(x)(isempty(x))||(isstruct(x)));
    p.addParameter('validateScenes', false, @islogical);
    p.addParameter('visualizeCustomConeFundamentals', false, @islogical);

    % Execute the parser
    p.parse(varargin{:});
    computeInputConeMosaicResponses = p.Results.computeInputConeMosaicResponses;
    computeInputConeMosaicResponsesBasedOnConeExcitations = p.Results.computeInputConeMosaicResponsesBasedOnConeExcitations;
    computeInputConeMosaicResponsesBasedOnPhotocurrents = p.Results.computeInputConeMosaicResponsesBasedOnPhotocurrents;
    inspectInputConeMosaicResponses = p.Results.inspectInputConeMosaicResponses;
    computeMRGCMosaicResponses = p.Results.computeMRGCMosaicResponses;
    mRGCNonLinearityParams = p.Results.mRGCNonLinearityParams;
    visualizeResponse = p.Results.visualizeResponse;
    visualizeStimulusSequence = p.Results.visualizeStimulusSequence;
    visualizeStimulusOnMosaic = p.Results.visualizeStimulusOnMosaic;
    debugInputConeMosaicPcurrentResponse = p.Results.debugInputConeMosaicPcurrentResponse;
    validateScenes = p.Results.validateScenes;
    visualizeCustomConeFundamentals = p.Results.visualizeCustomConeFundamentals;

    if (computeInputConeMosaicResponses)
        if (computeInputConeMosaicResponsesBasedOnConeExcitations)
            pCurrentTemporalResolutionSeconds = [];
            osBiophysicalModelWarmUpTimeSeconds = [];
            osBiophysicalModelTimeStep = [];
    
            computeInputConeMosaicTTF(theMRGCMosaic, theOI, TTFparamsStruct, ...
                theInputConeMosaicTTFResponsesFullFileName, ...
                pCurrentTemporalResolutionSeconds, osBiophysicalModelWarmUpTimeSeconds, osBiophysicalModelTimeStep, ...
                visualizeResponse, visualizeStimulusSequence, visualizeStimulusOnMosaic, debugInputConeMosaicPcurrentResponse, ...
                validateScenes, visualizeCustomConeFundamentals, ...
                'cone excitations');
        end
    
    
        if (computeInputConeMosaicResponsesBasedOnPhotocurrents) && ...
           (~isempty(mRGCNonLinearityParams)) && (strcmp(mRGCNonLinearityParams.type, 'photocurrent'))
    
            osBiophysicalModelWarmUpTimeSeconds = mRGCNonLinearityParams.osBiophysicalModelWarmUpTimeSeconds;
            osBiophysicalModelTimeStep = mRGCNonLinearityParams.osBiophysicalModelTemporalResolutionSeconds;
            pCurrentTemporalResolutionSeconds = mRGCNonLinearityParams.pCurrentTemporalResolutionSeconds;
    
            computeInputConeMosaicTTF(theMRGCMosaic, theOI, TTFparamsStruct, ...
                theInputConeMosaicTTFResponsesFullFileName, ...
                pCurrentTemporalResolutionSeconds, osBiophysicalModelWarmUpTimeSeconds, osBiophysicalModelTimeStep, ...
                visualizeResponse, visualizeStimulusSequence, visualizeStimulusOnMosaic, ...
                debugInputConeMosaicPcurrentResponse, inspectInputConeMosaicResponses, ...
                validateScenes, visualizeCustomConeFundamentals, ...
                'photocurrents');
        end
    end % computeInputConeMosaicResponses

    if (computeMRGCMosaicResponses)
        error('MRGCmosaic TTF computation not implemented yet')
        % Here we will need
    end

end


function computeInputConeMosaicTTF(theMRGCMosaic, theOI, TTFparamsStruct, ...
    theInputConeMosaicTTFResponsesFullFileName, ...
    pCurrentTemporalResolutionSeconds, photocurrentModelWarmUpTimeSeconds, osTimeStep, ...
    visualizeResponse, visualizeStimulusSequence, visualizeStimulusOnMosaic, ...
    debugInputConeMosaicPcurrentResponse, inspectInputConeMosaicResponses, ...
    validateScenes, visualizeCustomConeFundamentals, theSignal)

    % Compute the cone excitations
    if (strcmp(theSignal,'cone excitations'))
        fprintf('Input cone mosaic TTF responses will be saved in: \n%s\n', theInputConeMosaicTTFResponsesFullFileName);

        temporalFrequenciesExamined = TTFparamsStruct.tfSupport;

        % Compute cone contrasts for desired chromaticity
        coneContrasts = ...
            visualStimulusGenerator.coneContrastsFromChromaticity(TTFparamsStruct.chromaticity);

        totalContrast = TTFparamsStruct.contrast;


        stimParams = struct(...
            'displayType', TTFparamsStruct.displayType, ...
            'backgroundChromaticity', TTFparamsStruct.backgroundChromaticity, ...
            'backgroundLuminanceCdM2', TTFparamsStruct.backgroundLuminanceCdM2, ...
            'contrast', totalContrast, ...
            'coneContrasts', coneContrasts, ...
            'sizeDegs', TTFparamsStruct.sizeDegs, ...
            'positionDegs', TTFparamsStruct.positionDegs, ...
            'resolutionDegs', TTFparamsStruct.stimulusResolutionDegs, ...
            'stimulusPositionDegs', TTFparamsStruct.stimulusPositionDegs, ...
            'stimulusSizeDegs', TTFparamsStruct.stimulusSizeDegs, ...
            'stimulusShape', TTFparamsStruct.stimulusShape, ...
            'temporalPhaseIncrementDegs', TTFparamsStruct.temporalPhaseIncrementDegs, ...
            'coneMosaicModulationBasedResponse',  true ...
            );

        % Generate presentation display, 20% luminance headroom
        viewingDistanceMeters = 4;

        thePresentationDisplay = visualStimulusGenerator.presentationDisplay(...
            theMRGCMosaic.inputConeMosaic.wave, ...
            TTFparamsStruct.stimulusResolutionDegs, ...
            viewingDistanceMeters, ...
            'displayType', TTFparamsStruct.displayType, ...
            'bitDepth', 20, ...
            'meanLuminanceCdPerM2', TTFparamsStruct.backgroundLuminanceCdM2, ...
            'luminanceHeadroom', TTFparamsStruct.displayLuminanceHeadroomPercentage);


        if (TTFparamsStruct.coneFundamentalsOptimizedForStimPosition)
            % Compute custom cone fundamentals for the mosaic's eccentricity
            maxConesNumForAveraging = 3;
            customConeFundamentals = visualStimulusGenerator.coneFundamentalsForPositionWithinConeMosaic(...
                theMRGCMosaic.inputConeMosaic, theOI, stimParams.positionDegs, stimParams.sizeDegs, maxConesNumForAveraging);
        end

       
        for iTF = 1:numel(temporalFrequenciesExamined)

            % Get stim params
            stimParams.temporalFrequencyHz = temporalFrequenciesExamined(iTF);
            stimParams.durationSeconds =  1.0/stimParams.temporalFrequencyHz;

            % Feedback
            fprintf('Computing input cone mosaic TTF (cone excitations-based) for TF = %2.3f Hz\n', ...
                stimParams.temporalFrequencyHz);

            % Generate the spatial modulation patterns for all spatial phases of the drifting grating
            [theFlickeringGratingSpatialModulationPatterns, spatialSupportDegs, theCurrentTemporalPhasesDegs, ...
                theCurrentTemporalSupportSeconds, theCurrentTemporalRamp] = visualStimulusGenerator.flickeringGratingModulationPatterns(stimParams);


            if (TTFparamsStruct.coneFundamentalsOptimizedForStimPosition)
                if (iTF == 1)
                    % Generate scenes for the different frames of the flickring grating and for the null stimulus
                    [theFlickeringGratingFrameScenes, theNullStimulusScene, ~, theConeFundamentalsStruct] = visualStimulusGenerator.stimulusFramesScenes(...
                        thePresentationDisplay, stimParams, theFlickeringGratingSpatialModulationPatterns, ...
                        'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
                        'customConeFundamentals', customConeFundamentals, ...
                        'announceEmployedConeFundamentals', true, ...
                        'validateScenes', validateScenes, ...
                        'visualizeCustomConeFundamentals', visualizeCustomConeFundamentals);
                else
                    [theFlickeringGratingFrameScenes, theNullStimulusScene] = visualStimulusGenerator.stimulusFramesScenes(...
                        thePresentationDisplay, stimParams, theFlickeringGratingSpatialModulationPatterns, ...
                        'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
                        'withPreviouslyComputedConeFundamentalsStruct', theConeFundamentalsStruct, ...
                        'announceEmployedConeFundamentals', true, ...
                        'validateScenes', ~true);
                end

            else
                % Generate scenes for the different frames of the drifting grating and for the null stimulus
                [theFlickeringGratingFrameScenes, theNullStimulusScene] = visualStimulusGenerator.stimulusFramesScenes(...
                    thePresentationDisplay, stimParams, theFlickeringGratingSpatialModulationPatterns, ...
                    'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
                    'validateScenes', ~true);
            end

            % Compute input cone mosaic response to this temporal frequency

            % Place the stimulus at desired position
            if (isempty(stimParams.stimulusPositionDegs))
                stimulusPosition = 'mosaic-centered';
            elseif (numel(stimParams.stimulusPositionDegs) == 2)
                stimulusPosition = stimParams.stimulusPositionDegs;
            else
                stimParams.stimulusPositionDegs
                error('stimulusPositionDegs must be either empty of a 2 element vector')
            end

            % Set the mosaic integration time equal to the duration of one stimulus frame
            theMRGCMosaic.inputConeMosaic.integrationTime = theCurrentTemporalSupportSeconds(2) - theCurrentTemporalSupportSeconds(1);

            % Compute cone mosaic responses to each stimulus frame
            [theCurrentInputConeMosaicTTFresponse, theConeMosaicNullResponse] = ...
                RGCMosaicConstructor.helper.simulateExperiment.inputConeMosaicResponseToStimulusFrameSequence(...
                theMRGCMosaic, theOI, theNullStimulusScene, theFlickeringGratingFrameScenes, ...
                stimulusPosition, stimParams.coneMosaicModulationBasedResponse, ...
                'visualizeResponse', visualizeResponse, ...
                'visualizeStimulusSequence', visualizeStimulusSequence, ...
                'visualizeStimulusOnMosaic', visualizeStimulusOnMosaic, ...
                'thePresentationDisplayForVisualizingOpticalSceneOrImage', thePresentationDisplay, ...
                'stimulusInfoString', sprintf('TF:%2.2f Hz', stimParams.temporalFrequencyHz));

            if (iTF == 1)
                theInputConeMosaicTTFresponsesAllConditions = nan(...
                    numel(temporalFrequenciesExamined), ...
                    numel(theCurrentTemporalPhasesDegs), ...
                    theMRGCMosaic.inputConeMosaic.conesNum, ...
                    'single');

                temporalSupportSecondsAllConditions = nan(...
                    numel(temporalFrequenciesExamined), ...
                    numel(theCurrentTemporalSupportSeconds), ...
                    'single');

                temporalPhasesDegsAllConditions = nan(...
                    numel(temporalFrequenciesExamined), ...
                    numel(theCurrentTemporalPhasesDegs), ...
                    'single');
        
                temporalRampAllConditions = nan(...
                    numel(temporalFrequenciesExamined), ...
                    numel(theCurrentTemporalRamp), ...
                    'single');
            end % iTF == 1

            theInputConeMosaicTTFresponsesAllConditions(iTF, 1:size(theCurrentInputConeMosaicTTFresponse,1), 1:size(theCurrentInputConeMosaicTTFresponse,2)) = ...
                theCurrentInputConeMosaicTTFresponse;
            temporalSupportSecondsAllConditions(iTF, 1:numel(theCurrentTemporalSupportSeconds)) = single(theCurrentTemporalSupportSeconds);
            temporalPhasesDegsAllConditions(iTF, 1:numel(theCurrentTemporalPhasesDegs)) = single(theCurrentTemporalPhasesDegs);
            temporalRampAllConditions(iTF, 1:numel(theCurrentTemporalRamp)) = single(theCurrentTemporalRamp);
        end % iTF

        theInputConeMosaicTTFresponses = theInputConeMosaicTTFresponsesAllConditions;
        stimParams.temporalSupportSeconds = temporalSupportSecondsAllConditions;
        stimParams.temporalPhasesDegs = temporalPhasesDegsAllConditions;
        stimParams.temporalRamp = temporalRampAllConditions;

        % Save results
        fprintf('Saving computed input cone mosaic TTF cone excitation responses to %s.\n', theInputConeMosaicTTFResponsesFullFileName);
        stimParams.temporalFrequencyHz = temporalFrequenciesExamined;

        save(theInputConeMosaicTTFResponsesFullFileName, ...
            'theMRGCMosaic', ...
            'TTFparamsStruct', ...
            'stimParams', ...
            'theInputConeMosaicTTFresponses', ...
            'theConeMosaicNullResponse', ...
            '-v7.3');
    end % if (strcmp(theSignal,'cone excitations'))

    % Compute the photocurrents
    if (strcmp(theSignal, 'photocurrents'))
        load(theInputConeMosaicTTFResponsesFullFileName, ...
            'theMRGCMosaic', ...
            'TTFparamsStruct', ...
            'stimParams', ...
            'theInputConeMosaicTTFresponses', ...
            'theConeMosaicNullResponse');

        % If we loaded cone modulation responses, we will need to transform then to cone excitations.
        % Here we are reshaping the null response to allow for this transform
        if (stimParams.coneMosaicModulationBasedResponse)
            theConeMosaicNullResponse = squeeze(theConeMosaicNullResponse);
            theConeMosaicNullResponse = reshape(theConeMosaicNullResponse, [1 numel(theConeMosaicNullResponse)]);
        end

        % Read in the non-NaN cone-excitation data for each frequency
        % and compute the corresponding photocurrent response
        temporalFrequenciesExamined = stimParams.temporalFrequencyHz;

        for iTF = 1:numel(temporalFrequenciesExamined)

            % Retrieve the response data for this TF
            theTemporalSupportSecondsForThisTF = squeeze(stimParams.temporalSupportSeconds(iTF,:));
            theInputConeMosaicExcitationResponsesForThisTF = squeeze(theInputConeMosaicTTFresponses(iTF,:,:));
            nanIndices = find(isnan(theTemporalSupportSecondsForThisTF));
            if (isempty(nanIndices))
                theTimeBins = 1:numel(theTemporalSupportSecondsForThisTF);
            else
                theTimeBins = 1:(nanIndices-1);
            end
            theTemporalSupportSecondsForThisTF = theTemporalSupportSecondsForThisTF(theTimeBins);
            theInputConeMosaicExcitationResponsesForThisTF = theInputConeMosaicExcitationResponsesForThisTF(theTimeBins,:);

            % If we loaded cone modulation responses, transform then to cone excitations.
            if (stimParams.coneMosaicModulationBasedResponse)
                % Transform from modulations to excitations
                %Rmod = (Recx-Ro)/Ro -> (Rmod+1)*Ro = RExc
                theInputConeMosaicExcitationResponsesForThisTF = theInputConeMosaicExcitationResponsesForThisTF  + 1;
                theInputConeMosaicExcitationResponsesForThisTF = ...
                    bsxfun(@times, theInputConeMosaicExcitationResponsesForThisTF, theConeMosaicNullResponse);
            end

            % Compute the photocurrent response
            plotTitle = sprintf('tf: %2.3f Hz', stimParams.temporalFrequencyHz(iTF));

            fprintf('Computing input cone mosaic TTF (photocurrent-based) for TF = %2.3f Hz\n', ...
                stimParams.temporalFrequencyHz(iTF));

            stimulusPeriodDuration = 1/stimParams.temporalFrequencyHz(iTF); 
            nWarmUpPeriods = max([1 ceil(photocurrentModelWarmUpTimeSeconds/stimulusPeriodDuration)]);
            

            % Compute photocurrents here
            [theInputConeMosaicPhotocurrentTemporalSupportSecondsForThisTF, ...
             theInputConeMosaicPhotocurrentResponsesForThisTF, ...
             theInputConeMosaicBackgroundPhotocurrentsForThisTF] = RGCMosaicAnalyzer.compute.photocurrentsForOneStimulusPeriod(...
                    theMRGCMosaic.eccentricityDegs, ...
                    theTemporalSupportSecondsForThisTF(1:end-1), ...   % Get rid of last point which is a repeat of the first point
                    theInputConeMosaicExcitationResponsesForThisTF(1:end-1,:), ...       % Get rid of last point which is a repeat of the first point
                    nWarmUpPeriods, ...
                    pCurrentTemporalResolutionSeconds, osTimeStep, ...
                    theMRGCMosaic.inputConeMosaic.coneTypes, ...
                    debugInputConeMosaicPcurrentResponse, ...
                    plotTitle);

            timeBinIndices = 1:size(theInputConeMosaicPhotocurrentResponsesForThisTF,1);
            assert(numel(timeBinIndices) == numel(theInputConeMosaicPhotocurrentTemporalSupportSecondsForThisTF), 'mismatch in time dimension');
            assert(size(theInputConeMosaicPhotocurrentResponsesForThisTF,2) == theMRGCMosaic.inputConeMosaic.conesNum, 'error in cone index dimension');

            if (iTF == 1)
                theInputConeMosaicPhotocurrentTTFresponsesAllConditions = nan(...
                    numel(temporalFrequenciesExamined), ...
                    numel(timeBinIndices), ...
                    theMRGCMosaic.inputConeMosaic.conesNum, ...
                    'single');

                temporalSupportSecondsPhotocurrentAllConditions = nan(...
                    numel(temporalFrequenciesExamined), ...
                    numel(timeBinIndices), ...
                    'single');

                theInputConeMosaicBackgroundPhotocurrentsAllConditions = nan(...
                    numel(temporalFrequenciesExamined), ...
                    theMRGCMosaic.inputConeMosaic.conesNum, ...
                    'single');
            end % iTF == 1

            % Get the data for this TF into the all conditions array at the appropriate time bins

            theInputConeMosaicPhotocurrentTTFresponsesAllConditions(iTF, timeBinIndices, :) = theInputConeMosaicPhotocurrentResponsesForThisTF;
            temporalSupportSecondsPhotocurrentAllConditions(iTF, timeBinIndices) = theInputConeMosaicPhotocurrentTemporalSupportSecondsForThisTF;
            theInputConeMosaicBackgroundPhotocurrentsAllConditions(iTF,:) = theInputConeMosaicBackgroundPhotocurrentsForThisTF;

            if (inspectInputConeMosaicResponses)
                % Plot cone excitation + photocurrents
                hFig = figure(1000); clf;
                set(hFig, 'Position', [10 10 1000 1000]);
                set(hFig, 'Name', sprintf('mosaic cone excitations and photocurrents for TF = %2.3f Hz', temporalFrequenciesExamined(iTF)));

                ax = subplot(2,1,1);
                p1 = plot(ax, theTemporalSupportSecondsForThisTF*1e3, theInputConeMosaicExcitationResponsesForThisTF, 'r-');
                if (iTF == 1)
                    XLims = [theTemporalSupportSecondsForThisTF(1) theTemporalSupportSecondsForThisTF(end)]*1e3;
                end
                title(ax, 'cone excitations')
                xlabel(ax, 'time (msec)');
                ylabel(ax, 'cone excitation response');
                set(ax, 'XLim', XLims);

                ax = subplot(2,1,2);
                p2 = plot(ax, squeeze(temporalSupportSecondsPhotocurrentAllConditions(iTF, timeBinIndices))*1e3, ...
                              squeeze(theInputConeMosaicPhotocurrentTTFresponsesAllConditions(iTF, timeBinIndices, :)), 'b-');
                title(ax, 'photocurrents')
                set(ax, 'XLim', XLims);
                xlabel(ax, 'time (msec)');
                ylabel(ax, 'cone photocurrent response (pAmps)');

                drawnow
            end

        end % iTF

        % Append photocurrent data to theInputConeMosaicTTFResponsesFullFileName
        theInputConeMosaicPhotocurrentTemporalSupportSeconds = temporalSupportSecondsPhotocurrentAllConditions;
        theInputConeMosaicPhotocurrentTTFresponses = theInputConeMosaicPhotocurrentTTFresponsesAllConditions;
        theInputConeMosaicBackgroundPhotocurrents = theInputConeMosaicBackgroundPhotocurrentsAllConditions;

        save(theInputConeMosaicTTFResponsesFullFileName, ...
            'theInputConeMosaicPhotocurrentTemporalSupportSeconds', ...
            'theInputConeMosaicPhotocurrentTTFresponses', ...
            'theInputConeMosaicBackgroundPhotocurrents', ...
            '-append');

        fprintf('Photocurrent TTF responses were appeded to %s\n', theInputConeMosaicTTFResponsesFullFileName);

    end %if (strcmp(theSignal, 'photocurrents'))
end