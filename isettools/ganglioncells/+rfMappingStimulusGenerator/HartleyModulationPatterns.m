function [H, lIndices, mIndices] = HartleyModulationPatterns(omega, stimSizeDegs, pixelSizeDegs)

    pixelsNum  = round(stimSizeDegs / pixelSizeDegs);
    pixelIndex = 0:(pixelsNum-1);
    [X,Y] = meshgrid(pixelIndex/pixelsNum, pixelIndex/pixelsNum);

    H = zeros((2*omega+1)^2, size(X,1), size(X,2), 'single');
    lIndices = zeros((2*omega+1)^2,1);
    mIndices = lIndices;
    sIndex = 0;
    for mIndex = 0:(2*omega)
        for lIndex = 0:(2*omega)
            a = 2*pi*((lIndex-omega)*X + (mIndex-omega)*Y);
            f = sin(a)+cos(a);
            sIndex = sIndex + 1;
            H(sIndex,:,:) = single(f);
            lIndices(sIndex) = lIndex-omega;
            mIndices(sIndex) = mIndex-omega;
        end
    end
    H = H / max(H(:));
    H = single(H);
end