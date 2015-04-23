function [coneMosaicImage,coneMosaicRawImage] = generateConeMosaicImage(cfaPattern, coneSize, style)
% [coneMosaicImage,coneMosaicRawImage] = generateConeMosaicImage(cfaPattern, coneSize, style)
%
% Produce a viewable image of an isetbio cone mosaic.
% There are two available styles:
%   'standard' - this is the style used in Benson and Brainard (2014)
%   'williams_roorda' - this is the style used in Williams/Roorda papers
% 
% See also: t_generateConeMosaicImage
%
% 4/21/15  ncp  Wrote it.
% 4/22/15  dhb  Tweaks for isetbio style mosaics

    % Convert isetbio cfa to the representation Nicolas used when he wrote
    % this.
    sz = size(cfaPattern);
    coneMosaicRawImage = zeros([sz 3]);
    for ii = 1:3
        index = find(cfaPattern == ii+1);
        temp = zeros(sz);
        temp(index) = ii;
        coneMosaicRawImage(:,:,ii) = temp;
    end
    
    % Do the work
    subMosaicL = makeSingleConeMosaic(coneMosaicRawImage(:,:,1), coneSize, style);
    subMosaicM = makeSingleConeMosaic(coneMosaicRawImage(:,:,2), coneSize, style);
    subMosaicS = makeSingleConeMosaic(coneMosaicRawImage(:,:,3), coneSize, style);
    coneMosaicImage = makeThreeConeMosaic(subMosaicL, subMosaicM, subMosaicS, style);
end

function coneMosaic = makeThreeConeMosaic(mosaicL, mosaicM, mosaicS, style)
    coneMosaic = ones(size(mosaicL,1), size(mosaicL,2), 3);
    
    switch (style)
        case 'standard'
            % make green a bit dimmer
            mosaicM = mosaicM * 0.75;
        case 'williams_roorda'
            % boost the brightness a bit
            mosaicL = mosaicL / 0.75;
            mosaicM = mosaicM / 0.75;
            mosaicS = mosaicS / 0.75;
            mosaicL(mosaicL > 1) = 1;
            mosaicM(mosaicM > 1) = 1;
            mosaicS(mosaicS > 1) = 1;
    end
    
    coneMosaic(:,:,1) = mosaicL;
    coneMosaic(:,:,2) = mosaicM;
    coneMosaic(:,:,3) = mosaicS;
    
   switch (style)
        case 'standard'
        % Make all black areas between cones white
        indices = all(coneMosaic < 0.01,3);
        coneMosaic(repmat(indices,[1,1,3])) = 1;
   end
end

function singleConeMosaic = makeSingleConeMosaic(pixelMosaic, coneSize, style)

switch (style)
    case 'standard'
        
        kernel = makeCircularKernel(coneSize/2);
    case 'williams_roorda'
        kernel = makeGaussianKernel(coneSize/2);
end

singleConeMosaic = upsampleImage(pixelMosaic, coneSize+1);
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
    kernel(r < 1) = 1;
end

function upsampledImage = upsampleImage(originalImage, upsampleFactor)
    offset = (upsampleFactor-1)/2+1;
    [nRows,nCols] = size(originalImage);
    upsampledImage = zeros(upsampleFactor*nRows,upsampleFactor*nCols);
    upsampledImage(offset:upsampleFactor:end,offset:upsampleFactor:end) = originalImage;
end

