%
% visualStimulusGenerator.flickeringGratingModulationPatterns(stimParams)
%
function [H, spatialSupportDegs, temporalPhasesDegs, temporalSupportSeconds, temporalRamp] = ...
    flickeringGratingModulationPatterns(stimParams)


    spatialSupportDegs = visualStimulusGenerator.spatialSupport(...
        stimParams.sizeDegs, stimParams.resolutionDegs);
    [X,Y] = meshgrid(spatialSupportDegs,spatialSupportDegs);

    if (isfield(stimParams, 'temporalFrequencyHz') && (isfield(stimParams, 'durationSeconds')))
        oneCycleDurationSeconds = 1.0/stimParams.temporalFrequencyHz
        if (~isempty(stimParams.frameDurationMsec))
            framesPerCycle = round(oneCycleDurationSeconds*1e3 / stimParams.frameDurationMsec)
            stimParams.temporalPhaseIncrementDegs = 360/framesPerCycle;
        end
        framesNumPerCycle = 360/stimParams.temporalPhaseIncrementDegs;
        frameDurationSeconds = oneCycleDurationSeconds / framesNumPerCycle;
        if (~isempty(stimParams.frameDurationMsec))
            fprintf('\nDesired frame duration: %2.3f (msec), actual: %2.3f (msec)\n', stimParams.frameDurationMsec, frameDurationSeconds*1e3)
        end

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
            Router = max(stimParams.stimulusSizeDegs);
            theSpatialPattern(R > Router) = 0;
        case 'annulus'
            Rinner = min(stimParams.stimulusSizeDegs);
            Router = max(stimParams.stimulusSizeDegs);
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
