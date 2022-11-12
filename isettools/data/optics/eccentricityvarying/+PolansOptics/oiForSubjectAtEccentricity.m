function [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = oiForSubjectAtEccentricity(subjectID, whichEye, ecc, ...
    pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin)

    % Parse input
    p = inputParser;
    p.addRequired('subjectID', @(x)(isscalar(x)&&(x>=0)&&(x<=10)));
    p.addRequired('whichEye', @(x)(ischar(x)&&(ismember(x,{PolansOptics.constants.rightEye}))));  % allow only right eye data - the paper does not have left eye data
    p.addRequired('ecc', @(x)(isnumeric(x)&&(numel(x) == 2)));
    p.addRequired('pupilDiamMM', @(x)(isscalar(x)));
    p.addRequired('wavelengthsListToCompute', @(x)(isnumeric(x)));
    p.addRequired('micronsPerDegree', @(x)(isscalar(x)));
    p.addParameter('inFocusWavelength', 550, @isscalar);
    p.addParameter('wavefrontSpatialSamples', 801, @isscalar)
    p.addParameter('subtractCentralRefraction', true, @islogical);
    p.addParameter('zeroCenterPSF', false, @islogical);
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
    p.addParameter('noLCA', false, @islogical);
    p.addParameter('refractiveErrorDiopters', 0, @isnumeric);
    p.parse(subjectID, whichEye, ecc, pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin{:});

    % Check pupil diameter. Setting sujectID to 0 gets diffraction limited
    % optics so we skip check in that case.
    if (subjectID ~= 0)
        if (p.Results.pupilDiamMM < 1 || p.Results.pupilDiamMM > 4)
            error('Polans pupil diameter must be between 1 and 4 mm');
        end
    end
    
    inFocusWavelength = p.Results.inFocusWavelength;
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
    subtractCentralRefraction = p.Results.subtractCentralRefraction;
    zeroCenterPSF = p.Results.zeroCenterPSF;
    flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
    upsampleFactor = p.Results.upsampleFactor;
    noLCA = p.Results.noLCA;
    refractiveErrorDiopters = p.Results.refractiveErrorDiopters;
    refractiveErrorMicrons = wvfDefocusDioptersToMicrons(refractiveErrorDiopters, pupilDiamMM);
    
    % Obtain z-coeffs at desired eccentricity
    if (subjectID == 0)
        zCoeffs = zCoeffsForSubjectAtEcc(subjectID, ecc, subtractCentralRefraction, refractiveErrorMicrons);
        zCoeffs = 0*zCoeffs;
        measurementPupilDiameterMM = pupilDiamMM;
    else
        zCoeffs = zCoeffsForSubjectAtEcc(subjectID, ecc, subtractCentralRefraction, refractiveErrorMicrons);
        measurementPupilDiameterMM = PolansOptics.constants.measurementPupilDiamMM;
    end
    
    % Compute PSF and WVF from z-Coeffs for the desired pupil and wavelenghts
    [thePSF, ~, ~,~, psfSupportMinutesX, psfSupportMinutesY, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             measurementPupilDiameterMM, ...
             pupilDiamMM, inFocusWavelength, false, ...
             'doNotZeroCenterPSF', ~zeroCenterPSF, ...
             'micronsPerDegree', micronsPerDegree, ...
             'flipPSFUpsideDown', flipPSFUpsideDown, ...
             'upsampleFactor', upsampleFactor, ...
             'noLCA', noLCA, ...
             'name', sprintf('Polans subject %d, eccentricity: %2.1f,%2.1f degs', subjectID, ecc(1), ecc(2)));
    
    % Remove wavelength-dependent defocus if noLCA is set
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


function  interpolatedZcoeffs = zCoeffsForSubjectAtEcc(subjectID, ecc, subtractCentralRefraction, refractiveErrorMicrons)

    % Get original z-coeffs at all measured eccentricities
    if (subjectID == 0)
        [zMap, zCoeffIndices] = PolansOptics.constants.ZernikeCoefficientsMap(1);
    else
        [zMap, zCoeffIndices] = PolansOptics.constants.ZernikeCoefficientsMap(subjectID);
    end
    
    % Interpolate zMap at desired ecc
    [X,Y] = meshgrid(...
        PolansOptics.constants.measurementHorizontalEccentricities, ...
        PolansOptics.constants.measurementVerticalEccentricities);
    
    % index of coeffs at ecc = (0,0)
    indexOfZeroEcc = (X==0) & (Y==0);
    
    interpolatedZcoeffs = zeros(1, 30);
    zCoeffsNum = size(zMap,3);
    for zIndex = 1:zCoeffsNum
         % Retrieve the XY spatial map for this z-coeff
         z2Dmap = squeeze(zMap(:,:,zIndex));
         
         % The 4-th z-coeff is defocus. Subtract central defocus from all
         % spatial positions
         if ((zCoeffIndices(zIndex) == 4) && (subtractCentralRefraction))
             z2Dmap = z2Dmap - z2Dmap(indexOfZeroEcc);
         end
         
         % Add refractive error
         if (zCoeffIndices(zIndex) == 4)
             z2Dmap = z2Dmap + refractiveErrorMicrons;
         end
         
         % Interpolate the XY map at the desired eccentricity.
         interpolatedZcoeffs(zCoeffIndices(zIndex)+1) = interp2(X,Y,z2Dmap, ecc(1), ecc(2));
    end
    
    if (subjectID == 0)
        interpolatedZcoeffs = 0*interpolatedZcoeffs;
    end
    
end