function cmParams = cMosaicParams()
% Return a struct with the cMosaic default parameters
%
% If you control the cMosaic parameters at creation, the constructor
% returns your parameters as a second output argument
%
%   [cm, params] = cMosaic( ... your parameter ...) ;
%
% This makes it easy to re-create the cMosaic
%
%  cmDuplicate = cMosaic(params);
%
% See also
%   coneMosaicRectParams

cmParams.name = 'cone mosaic';
cmParams.wave =  400:10:700;
cmParams.pigment = cPhotoPigment();
cmParams.macular = Macular();
cmParams.coneData = [];

% In some cases we use eccentricityDegs.  BW likes position better,
% and he enabled using it as an alternative.
% cmParams.positionDegs =  [0 0];
cmParams.eccentricityDegs =  [0 0];

cmParams.sizeDegs =  [0.4 0.4];
cmParams.whichEye =  'right eye';
cmParams.computeMeshFromScratch =  false;
cmParams.customMinRFspacing =  [];
cmParams.customRFspacingFunction =  [];
cmParams.customDegsToMMsConversionFunction =  [];
cmParams.customMMsToDegsConversionFunction =  [];
cmParams.visualizeMeshConvergence =  false;
cmParams.exportMeshConvergenceHistory =  false;
cmParams.maxMeshIterations =  100;
cmParams.micronsPerDegree =  [];

% Eccentricity related
cmParams.eccVaryingConeAperture =  true;
cmParams.eccVaryingConeBlur =  false;
cmParams.eccVaryingOuterSegmentLength =  true;
cmParams.eccVaryingMacularPigmentDensity =  true;
cmParams.eccVaryingMacularPigmentDensityDynamic =  false;
cmParams.anchorAllEccVaryingParamsToTheirFovealValues =  false;
cmParams.rodIntrusionAdjustedConeAperture = 0.6;

cmParams.coneCouplingLambda =  [];
cmParams.coneApertureModifiers =  struct('smoothLocalVariations',true);
cmParams.coneDiameterToSpacingRatio =  1.0;
cmParams.coneDensities =  [0.6 0.3 0.1 0.0];
cmParams.tritanopicRadiusDegs =  0.15;
cmParams.integrationTime =  5/1000;
cmParams.opticalImagePositionDegs =  'mosaic-centered';

% Compute reltaed
cmParams.noiseFlag =  'random';  % Or frozen
cmParams.randomSeed = [];
cmParams.useParfor =  true;

end