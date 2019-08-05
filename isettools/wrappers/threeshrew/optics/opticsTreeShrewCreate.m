function optics = opticsTreeShrewCreate(varargin)
% Create an optics structure for the TreeShrew eye
%
% Syntax:
%   [optics, wvf] = OPTICSTREESHREWCREATE(varargin)
%

% Get the default tree-shrew optics params
defaultParams = opticsTreeShrewDefaultParams();

p = inputParser;
p.addParameter('name', '', @ischar);
p.addParameter('opticsType', defaultParams.opticsType, @ischar);
p.addParameter('inFocusPSFsigmaMicrons', defaultParams.inFocusPSFsigmaMicrons, @isnumeric);
p.addParameter('focalLengthMM', defaultParams.focalLengthMM, @isnumeric);
p.addParameter('pupilDiameterMM', defaultParams.pupilDiameterMM, @isnumeric);
p.addParameter('wavelengthSupport', 400:10:700, @isnumeric);
p.addParameter('maxSF', 20.0, @isnumeric);
p.addParameter('deltaSF', 0.1, @isnumeric);

% Parse input
p.parse(varargin{:});

% Set params
opticsName = p.Results.name;
inFocusPSFsigmaMicrons = p.Results.inFocusPSFsigmaMicrons;
pupilDiameterMM = p.Results.pupilDiameterMM;
wavelengthSupport = p.Results.wavelengthSupport;
opticsType = p.Results.opticsType;
deltaSF = p.Results.deltaSF;
maxSF = p.Results.maxSF;
focalLengthMM= p.Results.focalLengthMM;

% Compute microns per degree
focalLengthMeters = focalLengthMM / 1000;
posteriorNodalDistanceMM = focalLengthMM;
micronsPerDegree = posteriorNodalDistanceMM * 1000 * tand(1);
optics.micronsPerDegree = micronsPerDegree;

optics.type = 'optics';
optics.name = opticsName;
optics = opticsSet(optics, 'model', 'shiftInvariant');

switch lower(opticsType)
    case {'gaussian psf'}
        optics = opticsUpdateOTFUsingGaussianPSF(optics, inFocusPSFsigmaMicrons, ...
            maxSF, deltaSF, wavelengthSupport);
    otherwise
        error('Unknown optics type: ''%s''.', opticsType)
end % switch lower(opticsType)

% Set the focal length
optics = opticsSet(optics, 'focalLength', focalLengthMeters);

% Set the f-Number.
optics = opticsSet(optics, 'fnumber', focalLengthMeters*1000/pupilDiameterMM);

% no off-axis attenuation
optics = opticsSet(optics, 'offAxisMethod', 'skip');

% Pixel vignetting is off
optics.vignetting =  0;
end