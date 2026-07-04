function params = defaultParams(~)
% Return the default parameters for the cMosaic

% It is annoying that we have to update this function and the cMosaic
% constructor to keep them in sync.  Must be a better way.

params.name = 'cone mosaic';
params.macular = Macular();
params.pigment = cPhotoPigment();
params.wave = 400:10:700;

params.eccentricityDegs = [0,0];
params.sizeDegs = [0.4 0.4];
params.whichEye = 'left eye';

params.micronsPerDegree = [];
params.eccVaryingConeAperture = true;
params.eccVaryingOuterSegmentLength = true;
params.eccVaryingMacularPigmentDensity = true;

% New capabilities for more realistic but slower simulation
params.eccVaryingConeBlur = false;
params.eccVaryingMacularPigmentDensityDynamic  = false;

params.coneDensities = [0.6 0.3 0.1 0.0];
params.tritanopicRadiusDegs =  0.15;

params.noiseFlag = 'random';
params.randomSeed = [];
params.integrationTime = 5/1000;

% Parallel computations
params.useParfor = true;

end
