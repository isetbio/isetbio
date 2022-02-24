function cmParams = cMosaicParams()
% Return a struct with the cMosaic default parameters
%
% cMosaic itself when run with cMosaic('params',varargin)
% will return the parameters of the specifically created cmosaic.
% Shortly.
%
% See also
%   coneMosaicRectP

cmParams.name = 'cone mosaic';
cmParams.wave =  400:10:700;
cmParams.pigment = cPhotoPigment();
cmParams.macular = Macular();
cmParams.coneData = [];

% In some cases we use eccentricityDegs.  But BW likes this better.
cmParams.positionDegs =  [0 0];

cmParams.sizeDegs =  [0.4 0.4];
cmParams.whichEye =  'right eye';
cmParams.computeMeshFromScratch =  false;
cmParams.customMinRFspacing =  [];
cmParams.customRFspacingFunction =  [];
cmParams.customDegsToMMsConversionFunction =  [];
cmParams.customMMsToDegsConversionFunction =  [];
cmParams.visualizeMeshConvergence =  false;
cmParams.exportMeshConvergenceHistoryToFile =  false;
cmParams.maxMeshIterations =  100;
cmParams.micronsPerDegree =  [];
cmParams.eccVaryingConeAperture =  true;
cmParams.eccVaryingConeBlur =  false;
cmParams.eccVaryingOuterSegmentLength =  true;
cmParams.eccVaryingMacularPigmentDensity =  true;
cmParams.eccVaryingMacularPigmentDensityDynamic =  false;
cmParams.anchorAllEccVaryingParamsToTheirFovealValues =  false;
cmParams.coneCouplingLambda =  [];
cmParams.coneApertureModifiers =  struct('smoothLocalVariations',true);
cmParams.coneDiameterToSpacingRatio =  1.0;
cmParams.coneDensities =  [0.6 0.3 0.1 0.0];
cmParams.tritanopicRadiusDegs =  0.15;
cmParams.noiseFlag =  'random';  % Or frozen
cmParams.randomSeed = [];
cmParams.integrationTime =  5/1000;
cmParams.opticalImagePositionDegs =  'mosaic-centered';
cmParams.useParfor =  true;
end