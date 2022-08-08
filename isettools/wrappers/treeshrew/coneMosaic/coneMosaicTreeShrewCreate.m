function theConeMosaic = coneMosaicTreeShrewCreate(micronsPerDegree,varargin)
% Create a @coneMosaicHex object for the TreeShrew retina
%
% Syntax:
%   theConeMosaic = CONEMOSAICTREESHREWCREATE(varargin)
%
p = inputParser;
p.addParameter('fovDegs', [2 2], @isnumeric);
p.addParameter('spatialDensity', [0 0.5 0 0.5], @isnumeric);
p.addParameter('customLambda', 6.2, @isnumeric);
p.addParameter('customInnerSegmentDiameter', 6.0, @isnumeric);
p.addParameter('integrationTimeSeconds', 5/1000, @isnumeric);
p.addParameter('sConeMinDistanceFactor', 2.5, @isnumeric);
p.addParameter('resamplingFactor', 3, @isnumeric);
% Parse input
p.parse(varargin{:});

spatialDensity = p.Results.spatialDensity;
fovDegs = p.Results.fovDegs;
customLambda = p.Results.customLambda;
customInnerSegmentDiameter = p.Results.customInnerSegmentDiameter;
integrationTimeSeconds = p.Results.integrationTimeSeconds;
sConeMinDistanceFactor = p.Results.sConeMinDistanceFactor;
resamplingFactor = p.Results.resamplingFactor;

% Scale fovDegs: this is a hack because of the way the setSizeToFOV()
% is called in @coneMosaic, which assumes a 17mm focal length
fovDegs = fovDegs/((300/micronsPerDegree)^2);

if (spatialDensity(1) ~= 0)
    error('The first element in spatialDensity vector must be 0.');
end
if (spatialDensity(3) ~= 0)
    error('The third element in spatialDensity (M-cone density) vector must be 0.');
end

thePhotopigment = treeShrewPhotopigment();

theConeMosaic = coneMosaicHex(resamplingFactor, ...
    'fovDegs', fovDegs, ...
    'micronsPerDegree',micronsPerDegree, ...
    'integrationTime', integrationTimeSeconds, ...
    'pigment', thePhotopigment ,...
    'macular', macularPigmentTreeShrewCreate(thePhotopigment.wave), ...
    'customLambda', customLambda, ...
    'customInnerSegmentDiameter', customInnerSegmentDiameter, ...
    'spatialDensity', spatialDensity, ...
    'sConeMinDistanceFactor', sConeMinDistanceFactor, ...
    'sConeFreeRadiusMicrons', 0 ...
    );

end


