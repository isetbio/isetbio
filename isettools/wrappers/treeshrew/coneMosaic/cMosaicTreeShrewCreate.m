function theCMosaic = cMosaicTreeShrewCreate(varargin)
% Create a @cMosaic object for the TreeShrew retina
%
% Syntax:
%   theCMosaic = CMOSAICTREESHREWCREATE(varargin)
%
% Description:
%   This method first generates the old-style (@coneMosaicHex)
%   tree-shrew mosaic, it extracts the spatial cone data and 
%   generates the equivalent @cMosaic object
%
    p = inputParser;
    p.addParameter('fovDegs', [2 1], @isnumeric);
    p.addParameter('spatialDensity', [0 0.5 0 0.5], @isnumeric);
    p.addParameter('customLambda', 6.2, @isnumeric);
    p.addParameter('customInnerSegmentDiameter', 6.0, @isnumeric);
    p.addParameter('integrationTimeSeconds', 5/1000, @isnumeric);
    p.addParameter('sConeMinDistanceFactor', 2.5, @isnumeric);
    % Parse input
    p.parse(varargin{:});


    % Generate the treeshrew oi to retrieve the micronsPerDegree factor
    theOI = oiTreeShrewCreate('pupilDiameterMM', 2.0);
    micronsPerDegree = theOI.optics.micronsPerDegree;

    % Generate the old-style (@coneMosaicHex) cone mosaic
    theConeMosaicHex = coneMosaicTreeShrewCreate(...
        micronsPerDegree, ...
        'integrationTimeSeconds', p.Results.integrationTimeSeconds, ...
        'fovDegs', p.Results.fovDegs, ...
        'spatialDensity', p.Results.spatialDensity, ...
        'customLambda', p.Results.customLambda, ...
        'customInnerSegmentDiameter', p.Results.customInnerSegmentDiameter, ...
        'sConeMinDistanceFactor', p.Results.sConeMinDistanceFactor);
    

    % Generate a treeshrew @cPigment object for the tree tree shrew
    thePhotopigment = treeShrewCPhotopigment();

    % Generate a treeshrew @Macula object for the three shrew
    theMacularPigment = macularPigmentTreeShrewCreate(thePhotopigment.wave);

    % Generate a @cMosaic object for the tree shrew
    theCMosaic = cMosaic(...
        'coneData', theConeMosaicHex.coneData(), ...
        'micronsPerDegree', theConeMosaicHex.micronsPerDegree, ...
        'integrationTime', theConeMosaicHex.integrationTime, ...
        'pigment', thePhotopigment ,...
        'macular', theMacularPigment);
end