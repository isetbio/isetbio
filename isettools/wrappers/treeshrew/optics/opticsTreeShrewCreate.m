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
optics = opticsSet(optics, 'model', 'shiftInvariant');

switch lower(opticsType)
    case {'gaussian psf'}
        optics = opticsUpdateOTFUsingGaussianPSF(optics, inFocusPSFsigmaMicrons, ...
            maxSF, deltaSF, wavelengthSupport);
    case {'wvf'}
        optics = opticsFromTreeShrewZCoefs(pupilDiameterMM, wavelengthSupport, micronsPerDegree);
    otherwise
        error('Unknown optics type: ''%s''.', opticsType)
end % switch lower(opticsType)

% Set the optics name
optics.name = opticsName;

% Set the focal length
optics = opticsSet(optics, 'focalLength', focalLengthMeters);

% Set the f-Number.
optics = opticsSet(optics, 'fnumber', focalLengthMeters*1000/pupilDiameterMM);

% cos-4th off-axis attenuation
optics = opticsSet(optics, 'offAxisMethod', 'cos4th');

% Pixel vignetting is off
optics.vignetting =  0;
end

function optics = opticsFromTreeShrewZCoefs(pupilDiameterMM, wavelengthSupport, micronsPerDegree)

    % Reference: "Noninvasive imaging of the tree shrew eye: Wavefront
    % analysis and retinal imaging with correlative histology", Sajdak et
    % al 2019
    % 4 mm pupil used in measurements of Sajdak et al 2019
    measuredDiameterMM_TreeShrew = 4.0;
    
    % Defocusm with coefficient #4, (#5 in Matlab' indexing) is by far the dominant Zcoeffs in measurents of Sajdak et al 2019
    % So we are setting all the other coeffs to 0.
    % The data here are from Figure 2 of Sajdak et al 2019
    zCoeffs_TreeShrew = zeros(1,13);
    zCoeffs_TreeShrew(5) = -2.75;  % (microns)
    
    wvfP = wvfCreate(...
        'spatialsamples', 1001, ...
        'calc wavelengths', wavelengthSupport, ...
        'zcoeffs', zCoeffs_TreeShrew, ...
        'name', sprintf('treeshrew-%d', pupilDiameterMM), ...
        'umPerDegree', micronsPerDegree);
    
    wvfP = wvfSet(wvfP, 'measured pupil size', measuredDiameterMM_TreeShrew);
    wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMM);
    wvfP = wvfComputePSF(wvfP);
        
    optics = oiGet(wvf2oi(wvfP), 'optics');
end
