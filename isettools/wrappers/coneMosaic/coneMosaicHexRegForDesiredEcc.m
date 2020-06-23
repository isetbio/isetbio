% Create a regular hex cone mosaic with params for a desired eccentricity
%
% Description:
%    Create a regular hex cone mosaic with params for a desired
%    eccentricity. The mosaic is actually positioned at (0,0) but its
%    cones are spaced and have aperture appropriate for the desired
%    eccentricity.
% 
% See Also:
%   wrappers/scene/compute/sceneFromROI
%   ISETBIO LiveScripts repository: ls_computingWithEccentricityVaryingHexMosaics.mlx

% History:
%    7/9/19  NPC  ISETBIO Team, Copyright 2019

function theConeMosaic = coneMosaicHexRegForDesiredEcc(varargin)

% Obtain some default params
c = coneMosaicHex(1);
defaultMicronsPerDegree = c.micronsPerDegree;
defaultSconeFreeRadiusMicrons = c.sConeFreeRadiusMicrons;

% Params for a 7% s-cone population with a semirandom S-cone arrangement
spatialDensity = [0.5 0.25 0.15];
SconeMinDistanceFactor = 2;


p = inputParser;
p.addParameter('whichEye', 'left', @(x)ismember(x, {'left', 'right'}));   % which eye
p.addParameter('eccRadiusDegs', 30, @isnumeric);                          % radial eccentricity in degs
p.addParameter('eccAngleDegs', 45, @isnumeric);                           % angular eccentricity in degs
p.addParameter('fovDegs', [5 5], @isnumeric);                             % field of view in degs
p.addParameter('spatialDensity', spatialDensity, @isnumeric);             % density of L/M/S cones
p.addParameter('sConeMinDistanceFactor', SconeMinDistanceFactor, @isnumeric);
p.addParameter('integrationTimeSeconds', 5/1000, @isnumeric);
p.addParameter('resamplingFactor', 3, @isnumeric);

% Parse input
p.parse(varargin{:});
whichEye = p.Results.whichEye;
eccRadiusDegs = abs(p.Results.eccRadiusDegs);
eccAngleDegs = p.Results.eccAngleDegs;
fovDegs = p.Results.fovDegs;
integrationTimeSeconds = p.Results.integrationTimeSeconds;
resamplingFactor = p.Results.resamplingFactor;
sConeMinDistanceFactor = p.Results.sConeMinDistanceFactor;
spatialDensity = p.Results.spatialDensity;
spatialDensity = [0 spatialDensity(1) spatialDensity(2) spatialDensity(3)];


eccRadiusMeters = eccRadiusDegs*defaultMicronsPerDegree*1e-6;
[spacingMeters, apertureMeters, densityConesPerMM2] = ...
    coneSizeReadData('whichEye', whichEye, ...
    'eccentricity', eccRadiusMeters, ...
    'angle', eccAngleDegs);

coneSpacingMicrons = spacingMeters * 1e6;
coneApertureMicrons = apertureMeters * 1e6;

if (eccRadiusDegs > 0.0)
    sConeFreeRadiusMicrons = [];
else
    sConeFreeRadiusMicrons = defaultSconeFreeRadiusMicrons;
end

theConeMosaic = coneMosaicHex(resamplingFactor, ...
    'fovDegs', fovDegs, ...
    'integrationTime', integrationTimeSeconds, ...
    'customLambda', coneSpacingMicrons, ...
    'customInnerSegmentDiameter', coneApertureMicrons, ...
    'spatialDensity', spatialDensity, ...
    'sConeMinDistanceFactor', sConeMinDistanceFactor, ...
    'sConeFreeRadiusMicrons', sConeFreeRadiusMicrons ...
    );

end

