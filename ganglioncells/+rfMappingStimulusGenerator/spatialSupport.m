function x = spatialSupport(stimSize, pixelSize)
    assert(isscalar(stimSize),  'stimSize must be a scalar');
    assert(isscalar(pixelSize), 'pixelSize must be a scalar');

    pixelsNum  = round(stimSize / pixelSize);
    x = (1:pixelsNum)*pixelSize;
    x = x - mean(x);
 end
