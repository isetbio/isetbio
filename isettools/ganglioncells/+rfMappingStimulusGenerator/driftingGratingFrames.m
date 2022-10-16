function [H, spatialPhasesDegs] = driftingGratingFrames(stimParams)

    pixelsNum  = round(stimParams.stimSizeDegs / stimParams.pixelSizeDegs);
    x = (1:pixelsNum)*stimParams.pixelSizeDegs;
    x = x - mean(x);
    y = x;
    [X,Y] = meshgrid(x,y);

    spatialPhasesDegs = 0:stimParams.spatialPhaseIncrementDegs:360;
    framesNum = numel(spatialPhasesDegs);

    H = zeros(framesNum, size(X,1), size(X,2));
    fX = stimParams.spatialFrequencyCPD * cosd(stimParams.orientationDegs);
    fY = stimParams.spatialFrequencyCPD * sind(stimParams.orientationDegs);
    for frameIndex = 1:framesNum
        a = 2*pi*(fX*X +fY*Y) + spatialPhasesDegs(frameIndex)/180*pi;
        H(frameIndex,:,:) = cos(a);
    end

    figure(1);
    for frameIndex = 1:framesNum
        imagesc(x,y, squeeze(H(frameIndex,:,:)));
        axis 'image'
        drawnow
        pause(0.3);
    end

end