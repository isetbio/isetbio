function [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = oiForSubjectAtEccentricity(subjectID, whichEye, ecc, ...
    pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin)
% OIFORSUBJECTATECCENTRICITY Computes the Optical Image (OI) and Point Spread Function (PSF) for a specific subject's
% measured ocular aberration at a given eccentricity and pupil size, using the Artal et al. (2012) dataset.
%
% Syntax:
% [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
%    oiForSubjectAtEccentricity(subjectID, whichEye, ecc, pupilDiamMM, ...
%    wavelengthsListToCompute, micronsPerDegree, [varargin])
%
% Description:
% This function interpolates the measured Zernike coefficients from the Artal dataset
% (stored in ArtalOptics.constants.ZernikeCoefficientsMap) to estimate the
% ocular aberrations at a specified eccentricity. It then computes the Point Spread
% Function (PSF) and constructs an ISETBio Optical Image (OI) structure.
%
% Inputs:
%   subjectID (scalar, >= 0, <= 130):
%     The ID of the subject from the Artal dataset. Use 0 for an ideal, aberration-free eye.
%   whichEye (char, 'right eye' or 'left eye'):
%     The eye ('right eye' or 'left eye') for which to compute the optics.
%   ecc (numeric, scalar or 2-element vector):
%     The horizontal (and optional vertical) eccentricity (in degrees) at which to compute the aberrations.
%     If scalar, it is treated as horizontal eccentricity, and vertical is 0.
%   pupilDiamMM (scalar, >= 1, <= 4):
%     The diameter (in millimeters) of the simulated exit pupil.
%   wavelengthsListToCompute (numeric vector):
%     A vector of wavelengths (in nm) for which to compute the PSF and OI.
%   micronsPerDegree (scalar):
%     The conversion factor from visual degrees to microns on the retina, necessary for the OI structure.
%
% Optional Key/Value Pairs:
%   'inFocusWavelength' (scalar, default: 550):
%     The wavelength (in nm) that is assumed to be perfectly in focus (best focus).
%   'wavefrontSpatialSamples' (scalar, default: 801):
%     The number of spatial samples for the wavefront map (e.g., 801x801). Must be an odd number.
%   'subtractCentralRefraction' (logical, default: true):
%     If true, the defocus at ecc = 0 deg is subtracted from all Zernike coefficients at all eccentricities.
%   'zeroCenterPSF' (logical, default: false):
%     If true, the computed PSF is spatially centered around the peak.
%   'flipPSFUpsideDown' (logical, default: false):
%     If true, the PSF is flipped vertically (sometimes necessary for display/consistency).
%   'upsampleFactor' (scalar, default: []):
%     Optional upsampling factor for the PSF computation.
%   'noLCA' (logical, default: false):
%     If true, disables Longitudinal Chromatic Aberration (LCA).
%   'refractiveErrorDiopters' (numeric, default: 0):
%     Additional refractive error (in Diopters) to add as a defocus term.
%
% Outputs:
%   theOI (struct):
%     The ISETBio Optical Image structure, containing the calculated optics.
%   thePSF (matrix):
%     The Point Spread Function (PSF) as a 3D matrix (row x col x wavelength).
%   psfSupportMinutesX (vector):
%     The spatial support of the PSF along the X-dimension (in minutes of arc).
%   psfSupportMinutesY (vector):
%     The spatial support of the PSF along the Y-dimension (in minutes of arc).
%   psfSupportWavelength (vector):
%     A vector containing the wavelengths (in nm) for which the PSF was computed.
%   zCoeffs (vector):
%     The interpolated Zernike coefficients (microns) used to generate the optics (J_0 to J_30).
%
% See also: 
%   wvf2oiSpecial, computePSFandOTF, zCoeffsForSubjectAtEcc
%

%% Parse input
p = inputParser;
p.addRequired('subjectID', @(x)(isscalar(x)&&(x>=0)&&(x<=130)));
p.addRequired('whichEye', @(x)(ischar(x)&&(ismember(x,{ArtalOptics.constants.rightEye, ArtalOptics.constants.leftEye}))));
p.addRequired('ecc', @(x)(isnumeric(x)&&((isscalar(x))||(numel(x) == 2))));
p.addRequired('pupilDiamMM', @(x)(isscalar(x)&&(x>=1)&&(x<=4)));
p.addRequired('wavelengthsListToCompute', @(x)(isnumeric(x)));
p.addRequired('micronsPerDegree', @(x)(isscalar(x)));
p.addParameter('inFocusWavelength', 550, @isscalar);
p.addParameter('wavefrontSpatialSamples', 801, @isscalar)
p.addParameter('subtractCentralRefraction', true, @islogical);
p.addParameter('zeroCenterPSF', false, @islogical);
p.addParameter('flipPSFUpsideDown', false, @islogical);
p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(isscalar(x))&&(x>0))));
p.addParameter('noLCA', false, @islogical);
p.addParameter('refractiveErrorDiopters', 0, @isnumeric);
p.parse(subjectID, whichEye, ecc, pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin{:});

inFocusWavelength = p.Results.inFocusWavelength;
wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
subtractCentralRefraction = p.Results.subtractCentralRefraction;
zeroCenterPSF = p.Results.zeroCenterPSF;
upsampleFactor = p.Results.upsampleFactor;
flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
noLCA = p.Results.noLCA;

% Additional refractive error, if desired
refractiveErrorDiopters = p.Results.refractiveErrorDiopters;
refractiveErrorMicrons = wvfDefocusDioptersToMicrons(refractiveErrorDiopters, pupilDiamMM);

% Obtain z-coeffs at desired eccentricity
zCoeffs = zCoeffsForSubjectAtEcc(subjectID, whichEye, ecc(1), subtractCentralRefraction, refractiveErrorMicrons);

if (subjectID == 0)
    measurementPupilDiameterMM = pupilDiamMM;
else
    measurementPupilDiameterMM = ArtalOptics.constants.measurementPupilDiamMM;
end

if (isempty(zCoeffs))
    theOI = [];
    thePSF = [];
    psfSupportMinutesX = [];
    psfSupportMinutesY = [];
    psfSupportWavelength = [];
    return;
end

%% Compute PSF and WVF from z-Coeffs for the desired pupil and wavelenghts
[thePSF, ~, ~,~, psfSupportMinutesX, psfSupportMinutesY, theWVF] = ...
    computePSFandOTF(zCoeffs, ...
    wavelengthsListToCompute, wavefrontSpatialSamples, ...
    measurementPupilDiameterMM, ...
    pupilDiamMM, inFocusWavelength, false, ...
    'doNotZeroCenterPSF', ~zeroCenterPSF, ...
    'micronsPerDegree', micronsPerDegree, ...
    'flipPSFUpsideDown', flipPSFUpsideDown, ...
    'upsampleFactor',  upsampleFactor, ...
    'noLCA', noLCA, ...
    'name', sprintf('Artal subject %d, eccentricity: %2.1f,%2.1f degs', subjectID, ecc(1), 0));

% Remove wavelength-dependent defocus if noLCA is set.
%
% This also removes a wavelength dependent effect of diffraction,
% which I don't think is what we want to do.
% if (noLCA)
%     % Set all PSFs to the PSF at the in-focus wavelenth
%     [~,wTarget] = min(abs(wavelengthsListToCompute-inFocusWavelength));
%     targetPSF = thePSF(:,:,wTarget);
%     for waveIndex = 1:numel(wavelengthsListToCompute)
%         theWVF.psf{waveIndex} = targetPSF;
%         thePSF(:,:,waveIndex) = targetPSF;
%     end
% end

% Generate the OI from the wavefront map
theOI = wvf2oiSpecial(theWVF, micronsPerDegree, pupilDiamMM);

psfSupportWavelength = wavelengthsListToCompute;
end


function  interpolatedZcoeffs = zCoeffsForSubjectAtEcc(subjectID, whichEye, ecc, subtractCentralRefraction, refractiveErrorMicrons)

% Get original z-coeffs at all measured eccentricities
if (subjectID == 0)
    [zMap, zCoeffIndices] = ArtalOptics.constants.ZernikeCoefficientsMap(1, whichEye);
else
    [zMap, zCoeffIndices] = ArtalOptics.constants.ZernikeCoefficientsMap(subjectID, whichEye);
end

% Check that the coefficients do not contain NaNs
atLeastOneNaNCoefficient = any(isnan(zMap(:)));

% Check that the coefficients are not all 0
allZeroCoefficients = (sum(abs(zMap(:))) == 0);

if (atLeastOneNaNCoefficient || allZeroCoefficients)
    interpolatedZcoeffs = [];
    return;
end

zCoeffsNum = size(zMap,2);

% Index of coeffs at ecc = (0,0)
indexOfZeroEcc = find(ArtalOptics.constants.measurementHorizontalEccentricities==0);

interpolatedZcoeffs = zeros(1, 30);

% The Artal Z-coeff indices start from 0, not 1
theDefocusZcoeffIndex = wvfOSAIndexToVectorIndex('defocus') - 1;

for zIndex = 1:zCoeffsNum
    % Retrieve the X spatial map for this z-coeff
    z1Dmap = squeeze(zMap(:,zIndex));

    % Subtract central defocus from all spatial positions
    if ((zCoeffIndices(zIndex) == theDefocusZcoeffIndex) && (subtractCentralRefraction))
        z1Dmap = z1Dmap - z1Dmap(indexOfZeroEcc);
    end

    % Add refractive error
    if (zCoeffIndices(zIndex) == theDefocusZcoeffIndex)
        if (abs(refractiveErrorMicrons)>0)
            if (subtractCentralRefraction)
                fprintf('Original defocus at eccX = 0 (central-refraction subtracted) (microns): %f, added defocus (microns): %f\n', z1Dmap(indexOfZeroEcc), refractiveErrorMicrons);
            else
                fprintf('Original defocus at eccX = 0 (microns): %f, added defocus (microns): %f\n', z1Dmap(indexOfZeroEcc), refractiveErrorMicrons);
            end
        end
        z1Dmap = z1Dmap + refractiveErrorMicrons;
    end

    % Interpolate the XY map at the desired eccentricity.
    % Note the + 1 added to the z-coeff. This is to ensure that
    % we are addressing the correct indices expected by the wvf object
    interpolatedZcoeffs(zCoeffIndices(zIndex)+1) = interp1(...
        ArtalOptics.constants.measurementHorizontalEccentricities, z1Dmap, ecc(1));
end

if (subjectID == 0)
    % All zero coefficients
    interpolatedZcoeffs = 0*interpolatedZcoeffs;
    % Add refractive error
    interpolatedZcoeffs(wvfOSAIndexToVectorIndex('defocus')) = refractiveErrorMicrons;
end
end