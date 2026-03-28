%
% visualStimulusGenerator.flickeringGratingModulationPatterns(stimParams)
%
function [H, spatialSupportDegs, temporalPhasesDegs, temporalSupportSeconds, temporalRamp] = ...
    flickeringGratingModulationPatterns(stimParams)

    spatialSupportDegs = visualStimulusGenerator.spatialSupport(...
        stimParams.sizeDegs, stimParams.resolutionDegs);
    [X,Y] = meshgrid(spatialSupportDegs,spatialSupportDegs);


    if (isfield(stimParams, 'temporalFrequencyHz') && (isfield(stimParams, 'durationSeconds')))
        oneCycleDurationSeconds = 1.0/stimParams.temporalFrequencyHz;

        if (~isstruct(stimParams.temporalPhaseIncrementDegs))
            error('''stimParams.temporalPhaseIncrementDegs'' is expected to be a struct with 2 fields: ''val'' and ''maxFrameDurationMsec''.');
        end

        if (~isempty(stimParams.temporalPhaseIncrementDegs.maxFrameDurationMsec))
            framesNumPerCycle = 360/stimParams.temporalPhaseIncrementDegs.val;
            frameDurationMsec = 1e3 * oneCycleDurationSeconds / framesNumPerCycle;
            if (frameDurationMsec > stimParams.temporalPhaseIncrementDegs.maxFrameDurationMsec)
                fprintf('\n\nAt the examined TF (%2.2f Hz), with a phase increment of %2.1f degs, the frame duration (%2.2f msec) is > the specified maxFrameDurationMsec (%2.2f msec). Adjusting ...', ...
                    stimParams.temporalFrequencyHz, stimParams.temporalPhaseIncrementDegs.val, frameDurationMsec, stimParams.temporalPhaseIncrementDegs.maxFrameDurationMsec);
                
                desiredTemporalPhaseIncrementDegs = 360/(1e3*oneCycleDurationSeconds/stimParams.temporalPhaseIncrementDegs.maxFrameDurationMsec);
            else
                desiredTemporalPhaseIncrementDegs = 360/(1e3*oneCycleDurationSeconds/frameDurationMsec);
            end

            stimParams.temporalPhaseIncrementDegs = 360/round(360/desiredTemporalPhaseIncrementDegs);
            framesNumPerCycle = 360/stimParams.temporalPhaseIncrementDegs;
            frameDurationMsec = 1e3 * oneCycleDurationSeconds / framesNumPerCycle;
            correspondingTemporalFrequency = 1/(frameDurationMsec/1e3*framesNumPerCycle);

            fprintf('\nUsing a temporal phase increment of %2.2f degs, which results in %d frames/cycle and a frame duration of %2.1f msec with an actual TF of %2.2f Hz (specified TF: %2.2f Hz)\n', ...
                    stimParams.temporalPhaseIncrementDegs, framesNumPerCycle, frameDurationMsec, correspondingTemporalFrequency, stimParams.temporalFrequencyHz);
        else
            stimParams.temporalPhaseIncrementDeg = stimParams.temporalPhaseIncrementDeg.val;
            framesNumPerCycle = 360/stimParams.temporalPhaseIncrementDegs;
        end


        frameDurationSeconds = oneCycleDurationSeconds / framesNumPerCycle;
        framesNum = 1+round(stimParams.durationSeconds / frameDurationSeconds);
        temporalPhasesDegs = (0:(framesNum-1)) * stimParams.temporalPhaseIncrementDegs;
    else
        % Only generate one cycle
        temporalPhasesDegs = 0:stimParams.temporalPhaseIncrementDegs:360;
        framesNum = numel(temporalPhasesDegs);
        % Since we have no temporal information, arbitrarily assume the
        % stimulus is 1 second long
        frameDurationSeconds = 1.0/framesNum;
    end

    temporalPhasesDegs = mod(temporalPhasesDegs, 360);

    if (~(isfield(stimParams, 'contrastMask'))) || (isempty(stimParams.contrastMask))
        stimParams.contrastMask = X*0+1;
    end


    R = sqrt(X.^2 + Y.^2);
    theSpatialPattern = ones(size(X,1), size(X,2));

    switch (stimParams.stimulusShape)
        case 'spot'
            Router = 0.5*max(stimParams.stimulusSizeDegs);
            theSpatialPattern(R > Router) = 0;
        case 'annulus'
            Rinner = 0.5*min(stimParams.stimulusSizeDegs);
            Router = 0.5*max(stimParams.stimulusSizeDegs);
            if (Rinner == Router)
                error('Annulus stimulus must have 2 different sizes, inner, outer radii');
            end
            theSpatialPattern((R<Rinner)|(R>Router)) = 0;
        otherwise
            error('Unknown stimulusShape: ''%s''. Only know how to make ''spot'', and ''annulus''.', stimParams.stimulusShape)
    end % switch

    H = zeros(framesNum, size(X,1), size(X,2));
    temporalSupportSeconds = zeros(1, framesNum);

    for frameIndex = 1:framesNum
        a = temporalPhasesDegs(frameIndex)/180*pi;
        thePattern = sin(a) * theSpatialPattern;
        H(frameIndex,:,:) = thePattern .* stimParams.contrastMask;
        temporalSupportSeconds(frameIndex) = (frameIndex-1)*frameDurationSeconds;
    end

    if (isfield(stimParams, 'temporalEnvelopeTau'))
        tt = temporalSupportSeconds - mean(temporalSupportSeconds);
        temporalRamp = exp(-0.5 * (tt/ stimParams.temporalEnvelopeTau) .^ 4);
        for frameIndex = 1:framesNum
            H(frameIndex,:,:)  = H(frameIndex,:,:) * temporalRamp(frameIndex);
        end
    else
        temporalRamp = [1];
    end
