function theConeMosaic = coneMosaicTreeShrewCreate(micronsPerDegree,varargin)
% Create a coneMosaic object for the TreeShrew retina
%
% Syntax:
%   theConeMosaic = CONEMOSAICTREESHREWCREATE(varargin)
%
p = inputParser;
p.addParameter('fovDegs', [2 2], @isnumeric);
p.addParameter('lConeDensity', 0.5, @isnumeric);
p.addParameter('customLambda', 8, @isnumeric);
p.addParameter('customInnerSegmentDiameter', 7, @isnumeric);
p.addParameter('integrationTimeSeconds', 5/1000, @isnumeric);
p.addParameter('sConeMinDistanceFactor', 2, @isnumeric);
% Parse input
p.parse(varargin{:});

lConeDensity = p.Results.lConeDensity;
fovDegs = p.Results.fovDegs;
customLambda = p.Results.customLambda;
customInnerSegmentDiameter = p.Results.customInnerSegmentDiameter;
integrationTimeSeconds = p.Results.integrationTimeSeconds;
sConeMinDistanceFactor = p.Results.sConeMinDistanceFactor;

% Treeshrew-specific scaling
treeShrewScaling = 300/micronsPerDegree;


% Spatial densities of cones in the treeshrew
spatialDensity = [...
    0 ...
    lConeDensity ...   % L-cones
    0 ...              % M-cones
    1-lConeDensity ... %S-cones
];

thePhotopigment = treeShrewPhotopigment();

theConeMosaic = coneMosaicHex(7, ...
    'fovDegs', fovDegs/(treeShrewScaling^2), ...
    'micronsPerDegree',micronsPerDegree, ...
    'integrationTime', integrationTimeSeconds, ...
    'pigment', thePhotopigment ,...
    'macular', treeShrewMacularPigment(thePhotopigment.wave), ...
    'customLambda', customLambda, ...
    'customInnerSegmentDiameter', customInnerSegmentDiameter, ...
    'spatialDensity', spatialDensity, ...
    'sConeMinDistanceFactor', sConeMinDistanceFactor, ...
    'sConeFreeRadiusMicrons', 0 ...
    );
end

function theMacularPigment = treeShrewMacularPigment(wavelength)
% Generate threeshrew-specific macular pigment  (absent, so zero density)
theMacularPigment = Macular(...
    'wave', wavelength, ...
    'density', 0);
end

