function [H, spatialPhasesDegs, temporalSupportSeconds, temporalRamp] = driftingGratingFrames(stimParams)

    x = rfMappingStimulusGenerator.spatialSupport(...
        stimParams.stimSizeDegs, stimParams.pixelSizeDegs);
    y = x;
    [X,Y] = meshgrid(x,y);

    if (isfield(stimParams, 'temporalFrequencyHz') && (isfield(stimParams, 'durationSeconds')))
        oneCycleDurationSeconds = 1.0/stimParams.temporalFrequencyHz;
        framesNumPerCycle = 360/stimParams.spatialPhaseIncrementDegs;
        frameDurationSeconds = oneCycleDurationSeconds / framesNumPerCycle;
        framesNum = round(stimParams.durationSeconds / frameDurationSeconds);
        spatialPhasesDegs = (0:(framesNum-1)) * stimParams.spatialPhaseIncrementDegs;
    else
        % Only generate one cycle
        spatialPhasesDegs = 0:stimParams.spatialPhaseIncrementDegs:360;
        framesNum = numel(spatialPhasesDegs);
        % Since we have no temporal information, arbitrarily assume the
        % stimulus is 1 second long
        frameDurationSeconds = 1.0/framesNum;
    end

    spatialPhasesDegs = mod(spatialPhasesDegs, 360);

    if (~(isfield(stimParams, 'contrastMask'))) || (isempty(stimParams.contrastMask))
        stimParams.contrastMask = X*0+1;
    end

    H = zeros(framesNum, size(X,1), size(X,2));
    temporalSupportSeconds = zeros(1, framesNum);
    fX = stimParams.spatialFrequencyCPD * cosd(stimParams.orientationDegs);
    fY = stimParams.spatialFrequencyCPD * sind(stimParams.orientationDegs);


    for frameIndex = 1:framesNum
        a = 2*pi*(fX*X +fY*Y) + spatialPhasesDegs(frameIndex)/180*pi;
        thePattern = cos(a);
        H(frameIndex,:,:) = thePattern .* stimParams.contrastMask;
        temporalSupportSeconds(frameIndex) = (frameIndex-1)*frameDurationSeconds;
    end

    if (isfield(stimParams, 'temporalEnvelopeTau'))
        tt = temporalSupportSeconds - mean(temporalSupportSeconds);
        temporalRamp = exp(-0.5 * (tt/ stimParams.temporalEnvelopeTau) .^ 2);
        for frameIndex = 1:framesNum
            H(frameIndex,:,:)  = H(frameIndex,:,:)  * temporalRamp(frameIndex);
        end
    else
        temporalRamp = [1];
    end



end