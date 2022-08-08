function optics = opticsTreeShrewCreate(varargin)
% Create an optics structure for the TreeShrew eye
%
% Syntax:
%    [optics, wvf] = opticsTreeShrew;
%
% Description:
%    Set up an ISETBio tree shrew optics object (and wavefront optics
%    object).
%
% Inputs:
%   None.
%
% Outputs:
%   optics           - The optics object
%   wvf              - The wavefront object
%
% Optional key/value pairs.
%   'name'           - String.  Name to give the object. Default empty.
%   'opticsType'     - String. What type of optics
%                      'gaussian psf' - Use a Gaussian PSF
%                      'wvf' - PSF from Zernike coefficients
%   'whichShew'      - If wvf optics, which tree shrew data to use. 
%   'inFocusPSFsigmaMicrons' - If Gaussian PSF, sigma of PSF in microns.
%   'focalLengthMM'  - Focal length of tree shrew eye in mm.
%   'pupilDiameterMM' - Pupil diameter of tree shrew eye we're modeling in mm. Default 4 mm. 
%   'wavelengthSupport' - Wavelength support for the calculations
%   'measuredWavelength' - 840 nm light used in measurements of Sajdak et al 2019 (section 2.2).   
%                        However, we've set the defocus term in the table
%                        below to be best focus, which we are going to take
%                        as 550 nm.  So we set measuredWavelength to 550 by
%                        default.  LCA calculations are then done relative
%                        to this and should come out right.
%   'spatialSamples' - Linear number of pixels in representation of PSF.
%                      Default 801.
%   'psfSamplesPerMinute' - Linear pixel sampling density for the PSF.
%                      Default 0.05. This together with spatialSamples determines the size of
%                      the PSF support, and also the spatial sampling of the corresponding
%                      OTF.
%   'maxSF'          - Maximum spatial frequency.  Used in Gaussian PSF calculations
%                      to determine support for OTF calc that in the end
%                      drives everything. Default 20.
%   'deltaSF'        - Delta spatial frequency, also used in Gaussian PSF
%                      calcs. Default 0.1.
%
%   Default values for parameters not specified are set by
%   opticsTreeShrewDefaultParams.
% 
% See also: opticsTreeShrewDefaultParams.

% History:
%
%   08/06/22  dhb, eem  Finish up updating Zernikes for Sadjak et al.
%                       paper.

% Get the default tree-shrew optics params
defaultParams = opticsTreeShrewDefaultParams();

p = inputParser;
p.addParameter('name', '', @ischar);
p.addParameter('opticsType', defaultParams.opticsType, @ischar);
p.addParameter('whichShrew', defaultParams.whichShrew, @isnumeric);
p.addParameter('inFocusPSFsigmaMicrons', defaultParams.inFocusPSFsigmaMicrons, @isnumeric);
p.addParameter('focalLengthMM', defaultParams.focalLengthMM, @isnumeric);
p.addParameter('pupilDiameterMM', defaultParams.pupilDiameterMM, @isnumeric);
p.addParameter('wavelengthSupport', 400:10:700, @isnumeric);
p.addParameter('measuredWavelength',550, @isnumeric);
p.addParameter('spatialSamples', 801, @isnumeric);
p.addParameter('psfSamplesPerMinute', 0.05, @isnumeric);
p.addParameter('maxSF', 20.0, @isnumeric);
p.addParameter('deltaSF', 0.1, @isnumeric);

% Parse input
p.parse(varargin{:});

% Set params
opticsName = p.Results.name;
inFocusPSFsigmaMicrons = p.Results.inFocusPSFsigmaMicrons;
pupilDiameterMM = p.Results.pupilDiameterMM;
wavelengthSupport = p.Results.wavelengthSupport;
measuredWavelength = p.Results.measuredWavelength;
spatialSamples = p.Results.spatialSamples;
psfSamplesPerMinute = p.Results.psfSamplesPerMinute;
opticsType = p.Results.opticsType;
whichShrew = p.Results.whichShrew;
deltaSF = p.Results.deltaSF;
maxSF = p.Results.maxSF;
focalLengthMM = p.Results.focalLengthMM;

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
        optics = opticsFromTreeShrewZCoefs(whichShrew, pupilDiameterMM, wavelengthSupport, measuredWavelength, ...
            micronsPerDegree, spatialSamples, psfSamplesPerMinute);
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

function optics = opticsFromTreeShrewZCoefs(whichShrew, pupilDiameterMM, ...
    wavelengthSupport, measuredWavelength, micronsPerDegree, spatialSamples, ...
    psfSamplesPerMinute)

    % Reference: "Noninvasive imaging of the tree shrew eye: Wavefront
    % analysis and retinal imaging with correlative histology", Sajdak et
    % al 2019
    
    % 4 mm pupil used in measurements of Sajdak et al 2019
    measuredDiameterMM_TreeShrew = 4.0;
 
    % We have the tabulated ZCoeffs provided by Roorda from the paper
    % above.  Read those in here.
    %
    % Tag on a leading zero because the file does not contain the first
    % (piston) Zernike coefficient, but ISETBio expects it.
    zCoeffsFile = fullfile(isetbioRootPath,'data','datafiles','treeshrew','SajdakEtAlTreeShrewZernikes.xlsx');
    zCoeffsAll = cell2mat(readcell(zCoeffsFile,'Range','B2:L66'));
    zCoeffs_TreeShrew = [0 ; zCoeffsAll(:,whichShrew)];
   
    % Create the wavefront object
    wvfP = wvfCreate(...
        'spatialsamples', spatialSamples, ...
        'measured wl', measuredWavelength, ... 
        'calc wavelengths', wavelengthSupport, ...
        'zcoeffs', zCoeffs_TreeShrew, ...
        'name', sprintf('treeshrew-%d', pupilDiameterMM), ...
        'umPerDegree', micronsPerDegree, ...
        'customLCA', @treeShrewLCA);

    % Check pupil diameter OK and then set
    if (pupilDiameterMM > measuredDiameterMM_TreeShrew)
        error('Requested pupil diameter greated than measured diameter');
    end
    wvfP = wvfSet(wvfP, 'measured pupil size', measuredDiameterMM_TreeShrew);
    wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMM);
    wvfP = wvfSet(wvfP, 'ref psf sample interval', psfSamplesPerMinute);
    % Compute the PSF using the wavefront code
    wvfP = wvfComputePSF(wvfP);
    
    % Create the corresponding optics object
    optics = oiGet(wvf2oi(wvfP), 'optics');

    % Set the microns per degree in the optics as well
    optics.micronsPerDegree = micronsPerDegree;
end

function lcaDiopters = treeShrewLCAOLD(measuredWavelength, thisWavelength)
    % We dont have a model LCA for tree shrew yet.
    % Here we model is as the human LCA x 5
    % This creates an LCA difference between 840nm and 550nm of -4.5D
    % as per Sadjak et al (2019), "Noninvasive imaging of the tree shrew eye:
    % Wavefront analysis and retinal imaging with correlative histology"
    
    constant = 1.8859 - (0.63346 ./ (0.001 .* measuredWavelength - 0.2141));
    lcaDiopters = 1.8859 - constant - (0.63346 ./ (0.001 * thisWavelength - 0.2141));
    lcaDiopters = 5.15 * lcaDiopters;
end

