function optics = opticsTreeShrewCreate(varargin)
% Create an optics structure for the TreeShrew eye
%
% Syntax:
%   [optics, wvf] = OPTICSTREESHREWCREATE(varargin)
%
% Default params:
% The peak spatial frequency sensitivity is around 0.5 cyc/deg, 
% with a high-frequency cutoff of 2 c/deg (see "Spatial contrast sensitivity of the tree shrew". Petry HM, Fox R, Casagrande VA.)
% 
% From "Normal development of refractive state and ocular component dimensions 
% in the tree shrew (Tupaia belangeri)", by Thomas T. Norton, Neville A. McBrien, 1992
% anteriorFocalLengthMM = 4.35;
% posteriorNodalDistanceMM = posteriorFocalPoint - posteriorNodalPoint = 7.84-3.49 = 4.35;
%

p = inputParser;
p.addParameter('name', '', @ischar);
p.addParameter('opticsType', 'gaussian psf', @ischar);
p.addParameter('inFocusPSFsigmaMicrons', 20, @isnumeric);
p.addParameter('pupilDiameterMM', 4, @isnumeric);
p.addParameter('wavelengthSupport', 400:10:700, @isnumeric);
p.addParameter('maxSF', 10.0, @isnumeric);
p.addParameter('deltaSF', 0.02, @isnumeric);
p.addParameter('focalLengthMM', 4.35, @isnumeric);

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
        optics = opticsWithGaussianPSF(optics, inFocusPSFsigmaMicrons, ...
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


function optics = opticsWithGaussianPSF(optics, psfSigmaMicrons, maxSF, deltaSF, wavelengthSupport)
    % Generate spatial frequency support
    N = maxSF/deltaSF;
    fList = unitFrequencyList(N);  % Should we add fList(end+1) = -fList(1) ?? (NPC)
    fList = fList * maxSF;
    [xSfGridCyclesDeg,ySfGridCyclesDeg] = meshgrid(fList, fList);
    fSupportCyclesPerDeg(:, :, 1) = xSfGridCyclesDeg;
    fSupportCyclesPerDeg(:, :, 2) = ySfGridCyclesDeg;
    
    % Compute the sigma of the PSF in minutes of arc, 1 deg = 60 minutes or arc
    psfSigmaMinutes = psfSigmaMicrons / optics.micronsPerDegree * 60;
    
    % Generate the 2D OTF (single wavelength) based on a Gaussian PSF
    theOTF = otfFromGaussianPSF(psfSigmaMinutes, fSupportCyclesPerDeg);
    
    % Reshape OTF in isetbio format
    theOTFIsetbioFormat = ifftshift(theOTF);
    otfData = zeros(size(theOTFIsetbioFormat,1), size(theOTFIsetbioFormat,2), numel(wavelengthSupport));
    for iW = 1:numel(wavelengthSupport)
        otfData(:,:,iW) = theOTFIsetbioFormat;
    end
    
    % Compute OTF frequency support in cycles/mm
    mmPerDegree = optics.micronsPerDegree * 1e-3;
    fSupportCyclesPerMM = fSupportCyclesPerDeg / mmPerDegree;
    fxCyclesPerMM = fSupportCyclesPerMM(1, :, 1);
    fyCyclesPerMM = fSupportCyclesPerMM(:, 1, 2);
    
    optics = opticsSet(optics, 'otf data',otfData);
    optics = opticsSet(optics, 'otffx', fxCyclesPerMM(:)');
    optics = opticsSet(optics, 'otffy', fyCyclesPerMM(:)');
    optics = opticsSet(optics, 'otfWave', wavelengthSupport);
end

function theOTF = otfFromGaussianPSF(psfSigmaMinutes, fSupportCyclesPerDeg)
    xSfGridCyclesDeg = fSupportCyclesPerDeg(:, :, 1);
    ySfGridCyclesDeg = fSupportCyclesPerDeg(:, :, 2);
    [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(xSfGridCyclesDeg,ySfGridCyclesDeg);
    thePSF = gaussianPSF(xGridMinutes, yGridMinutes, psfSigmaMinutes);
    [~,~,theOTF] = PsfToOtf(xGridMinutes,yGridMinutes,thePSF);
end

function thePSF = gaussianPSF(xGridMinutes, yGridMinutes, sigmaMinutes)
    radius = sqrt(xGridMinutes.^2 + yGridMinutes.^2);
    thePSF = normpdf(radius,0,sigmaMinutes);
    % Unit volume
    thePSF = thePSF/sum(thePSF(:));
end