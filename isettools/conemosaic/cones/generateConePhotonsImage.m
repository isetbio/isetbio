function [photonsImageNormalized,photonsImageRaw] = generateConePhotonsImage(photons, coneSize, style)
% [photonImageNormalized,photonsImageRaw] = generateConePhotonsImage(photons, coneSize, style)
%
% Produce a viewable image of the photon response array.  These are meant
% to be paired with cone mosaic images produced by generateConeMosaicImage.
% There are two available styles:
%   'standard' - this is the style used in Benson and Brainard (2014)
%   'williams_roorda' - this is the style used in Williams/Roorda papers
% 
% See also: generateConePhotonsImage, t_generateConeMosaicImage
%
% 4/27/15  dhb  Wrote from NCP's generateConePhotonsImage.
    
    % Do the work
    factor = max(photons(:));
    photonsImageNormalized = makeSingleConeMosaic(photons/factor, coneSize, style);
    photonsImageRaw = factor*photonsImageNormalized;
    
end

function singleConeMosaic = makeSingleConeMosaic(pixelMosaic, coneSize, style)

switch (style)
    case 'standard'
        kernel = makeCircularKernel(coneSize/2);
    case 'williams_roorda'
        kernel = makeGaussianKernel(coneSize/2);
end

singleConeMosaic = upsampleImage(pixelMosaic, coneSize);
singleConeMosaic = conv2(singleConeMosaic, kernel, 'same');
end


function kernel = makeGaussianKernel(kernelSize)
    [X,Y] = meshgrid(-kernelSize:kernelSize,-kernelSize:kernelSize);
    sigma = kernelSize/2.2;
    kernel = exp(-0.5*(X.^2 + Y.^2)/sigma^2);
end

function kernel = makeCircularKernel(kernelSize)
    [X,Y] = meshgrid(-kernelSize:kernelSize,-kernelSize:kernelSize);
    kernel = zeros(size(X));
    r = sqrt(X.^2 + Y.^2)/kernelSize;
    kernel(r < 0.9) = 1;
end

function upsampledImage = upsampleImage(originalImage, upsampleFactor)
    offset = round(upsampleFactor/2+1);
    [nRows,nCols] = size(originalImage);
    upsampledImage = zeros(upsampleFactor*nRows+1,upsampleFactor*nCols+1);
    upsampledImage(offset:upsampleFactor:end,offset:upsampleFactor:end) = originalImage;
end



