function theConeMosaic = coneMosaicTreeShrewCreate(micronsPerDegree,varargin)
% Create a coneMosaic object for the TreeShrew retina
%
% Syntax:
%   theConeMosaic = CONEMOSAICTREESHREWCREATE(varargin)
%
p = inputParser;
p.addParameter('fovDegs', [0.5 0.5], @isnumeric);
p.addParameter('lConeDensity', 0.7, @isnumeric);
p.addParameter('coneSeparationMicrons', 6, @isnumeric);
p.addParameter('innerSegmentDiameterMicrons', 6, @isnumeric);
p.addParameter('integrationTimeSeconds', 50/1000, @isnumeric);
% Parse input
p.parse(varargin{:});
fovDegs = p.Results.fovDegs;
lConeDensity = p.Results.lConeDensity;
coneSeparationMicrons = p.Results.coneSeparationMicrons;
innerSegmentDiameterMicrons = p.Results.innerSegmentDiameterMicrons;
integrationTimeSeconds = p.Results.integrationTimeSeconds;

% Load treeshrew-specific cone photopigments and macular pigment
theConePhotopigments = treeShrewConePhotopigments();
theMacularPigment = treeShrewMacularPigment(theConePhotopigments.wave);

% Spatial densities of cones in the treeshrew
spatialDensity = [...
    0 ...
    lConeDensity ...   % L-cones
    0 ...              % M-cones
    1-lConeDensity ... %S-cones
];

theConeMosaic = coneMosaicHex(1, ...
    'micronsPerDegree', micronsPerDegree, ...
    'fovDegs', fovDegs, ...
    'integrationTime', integrationTimeSeconds, ...
    'pigment', theConePhotopigments ,...
    'macular', theMacularPigment, ...
    'customLambda', coneSeparationMicrons, ...
    'customInnerSegmentDiameter', innerSegmentDiameterMicrons, ...
    'spatialDensity', spatialDensity, ...
    'sConeMinDistanceFactor', 0, ...
    'sConeFreeRadiusMicrons', 0 ...
    );
end

function theMacularPigment = treeShrewMacularPigment(wavelength)
% Generate threeshrew-specific macular pigment  (absent, so zero density)
theMacularPigment = Macular(...
    'wave', wavelength, ...
    'density', 0);
end

function thePhotopigment = treeShrewConePhotopigments() 
% Load the threeshrew L and S-cone absorbance spectra
load('treeshrewConeAbsorbanceSpectra.mat', 'wavelength', 'data');

% Max values of optical densities
peakOpticalDensitiesLS = max(data,[],1);

normalizedAbsorbanceSpectra = zeros(numel(wavelength),3);
% L-cone normalized absorbance spectra
normalizedAbsorbanceSpectra(:,1) = data(:,1)/peakOpticalDensitiesLS(1);
% S-cone normalized absorbance spectra
normalizedAbsorbanceSpectra(:,3) = data(:,2)/peakOpticalDensitiesLS(2);

% Generate treeshrew-specific cone photopigments
thePhotopigment = photoPigment(...
    'wave', wavelength, ...
    'absorbance', normalizedAbsorbanceSpectra, ...
    'opticalDensity', [peakOpticalDensitiesLS(1) 0 peakOpticalDensitiesLS(2)], ...
    'peakEfficiency', 0.5*[1 1 1]);
end

