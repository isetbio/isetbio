function [H, lIndices, mIndices] = HartleyModulationPatterns(X,Y,omega)

    H = zeros((2*omega+1)^2, size(X,1), size(X,2));
    lIndices = zeros((2*omega+1)^2,1);
    mIndices = lIndices;
    sIndex = 0;
    for mIndex = 0:(2*omega)
        for lIndex = 0:(2*omega)
            a = 2*pi*((lIndex-omega)*X + (mIndex-omega)*Y);
            f = sin(a)+cos(a);
            sIndex = sIndex + 1;
            H(sIndex,:,:) = f;
            lIndices(sIndex) = lIndex-omega;
            mIndices(sIndex) = mIndex-omega;
        end
    end
    H = H / max(H(:));
end