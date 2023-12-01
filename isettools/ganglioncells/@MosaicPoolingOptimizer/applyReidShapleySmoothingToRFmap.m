function [smoothedRF, smoothingKernel, rfPixelSizeSamples] = applyReidShapleySmoothingToRFmap(spatialSupportDegsX, theRFmap, rfPixelsAcross)
    
    dx = (spatialSupportDegsX(2)-spatialSupportDegsX(1));
    rfPixelSizeDegs = (spatialSupportDegsX(end)-spatialSupportDegsX(1)+dx)/rfPixelsAcross;
    rfPixelSizeSamples = rfPixelSizeDegs/dx;
        
    xx = -rfPixelSizeSamples:rfPixelSizeSamples;
    [X, Y] = meshgrid(xx*dx,xx*dx);

    % Drop to 10%, at half pixel (ReidShapley, J. Neuroscience, 2002)
    sigma = 0.23*rfPixelSizeDegs;

    h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
    smoothingKernel = h / sum(h(:));

    smoothedRF = conv2(theRFmap, smoothingKernel, 'same');
end