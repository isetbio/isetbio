function cmParams = coneMosaicRectP()
% Return a struct with default coneMosaicRect parameters
%
%
% See also
%   coneMosaic
%

cmParams.name    = 'rectMosaic';
cmParams.pigment = photoPigment();
cmParams.macular = Macular();
cmParams.os      = [];
cmParams.center  = [0 0];
cmParams.wave    = (400:10:700);
cmParams.pattern = [];
cmParams.spatialDensity = [0 0.6 0.3 0.1];
cmParams.size    = [72 88];
cmParams.integrationTime  = 0.005;
cmParams.micronsPerDegree = 300;
cmParams.emPositions      = [0 0];
cmParams.apertureBlur     = false;
cmParams.noiseFlag        = 'random';
cmParams.eccentricityunits = 'm';
cmParams.whichEye = 'left';
cmParams.apertureBlur = false;
cmParams.useParFor = false;

end
