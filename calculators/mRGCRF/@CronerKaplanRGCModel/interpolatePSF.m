function [thePSFHR, thePSFsupportDegsHR] = interpolatePSF(thePSF, thePSFsupportDegs, upsampleFactor, paddingMarginDegs)
    
    % Upsample the PSF via linear interpolation
    thePSFsupportDegsHR = linspace(thePSFsupportDegs(1), thePSFsupportDegs(end), numel(thePSFsupportDegs)*upsampleFactor);
    [X2,Y2] = meshgrid(thePSFsupportDegsHR,thePSFsupportDegsHR);
    [X,Y] = meshgrid(thePSFsupportDegs, thePSFsupportDegs);
    thePSFHR = interp2(X,Y, thePSF, X2, Y2);
    thePSFHR = thePSFHR / sum(thePSFHR(:));
    
    % zero pad to minimize edge artifacts
    sampleSizeDegs = thePSFsupportDegsHR(2)-thePSFsupportDegsHR(1);
    paddingPixelsNum = round(paddingMarginDegs/sampleSizeDegs);
    minPSFsupportDegs = thePSFsupportDegs(1)-paddingPixelsNum*sampleSizeDegs;
    maxPSFsupportDegs = thePSFsupportDegs(end)+paddingPixelsNum*sampleSizeDegs;
    newPixelsNum = length(thePSFsupportDegsHR) + 2*paddingPixelsNum;
    
    % New PSF support
    thePSFsupportDegsHR = linspace(minPSFsupportDegs, maxPSFsupportDegs, newPixelsNum);
    
    % Zero pad
    tmp = thePSFHR;
    thePSFHR = zeros(newPixelsNum,newPixelsNum);
    thePSFHR(paddingPixelsNum+(1:size(tmp,1)), paddingPixelsNum+(1:size(tmp,2))) = tmp;
end