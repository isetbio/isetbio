function [H, spatialPhasesDegs] = driftingGratingFrames(stimParams)

    x = rfMappingStimulusGenerator.spatialSupport(...
        stimParams.stimSizeDegs, stimParams.pixelSizeDegs);
    y = x;
    [X,Y] = meshgrid(x,y);

    spatialPhasesDegs = 0:stimParams.spatialPhaseIncrementDegs:360;
    framesNum = numel(spatialPhasesDegs);

    if (~(isfield(stimParams, 'contrastMask'))) || (isempty(stimParams.contrastMask))
        stimParams.contrastMask = X*0+1;
    end

    H = zeros(framesNum, size(X,1), size(X,2));
    fX = stimParams.spatialFrequencyCPD * cosd(stimParams.orientationDegs);
    fY = stimParams.spatialFrequencyCPD * sind(stimParams.orientationDegs);
    for frameIndex = 1:framesNum
        a = 2*pi*(fX*X +fY*Y) + spatialPhasesDegs(frameIndex)/180*pi;
        thePattern = cos(a);
        H(frameIndex,:,:) = thePattern .* stimParams.contrastMask;
    end

end