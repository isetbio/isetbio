function [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = ...
    oiForSubjectAtEccentricity(subjectID, whichEye, ...
    pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin)

    % Parse input
    p = inputParser;
    p.addRequired('subjectID', @(x)(isscalar(x)&&(x>=0)&&(x<=130)));
    p.addRequired('whichEye', @(x)(ischar(x)&&(ismember(x,{ThibosOptics.constants.rightEye, ThibosOptics.constants.leftEye}))));
    p.addRequired('pupilDiamMM', @(x)(isscalar(x)&&(x>=1)&&(x<=max(ThibosOptics.constants.availableMeasurementPupilDiamsMM))));
    p.addRequired('wavelengthsListToCompute', @(x)(isnumeric(x)));
    p.addRequired('micronsPerDegree', @(x)(isscalar(x)));
    p.addParameter('measurementPupilDiameterMM', ThibosOptics.constants.measurementPupilDiamMM,  @(x)(isscalar(x) && ismember(x, ThibosOptics.constants.availableMeasurementPupilDiamsMM)));
    p.addParameter('inFocusWavelength', 550, @isscalar);
    p.addParameter('wavefrontSpatialSamples', 801, @isscalar)
    p.addParameter('subtractCentralRefraction', true, @islogical);
    p.addParameter('zeroCenterPSF', false, @islogical);
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
    p.addParameter('noLCA', false, @islogical);
    p.addParameter('refractiveErrorDiopters', 0, @isnumeric);
    p.parse(subjectID, whichEye, pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin{:});
  

    measurementPupilDiameterMM = p.Results.measurementPupilDiameterMM;
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
   
    % Obtain z-coeffs
    zCoeffs = zCoeffsForSubject(subjectID, whichEye, measurementPupilDiameterMM, subtractCentralRefraction, refractiveErrorMicrons);

    if (subjectID == 0)
        % Subject ID is a special case: diffraction-limited optics  possibly with some
        % refractive error (see zCoeffsForSubject)
        measurementPupilDiameterMM = pupilDiamMM;
    end

    if (isempty(zCoeffs))
        theOI = [];
        thePSF = []; 
        psfSupportMinutesX = []; 
        psfSupportMinutesY = [];
        psfSupportWavelength = [];
        return;
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
             'upsampleFactor',  upsampleFactor, ...
             'noLCA', noLCA, ...
             'name', sprintf('Thibos subject %d, eccentricity: %2.1f,%2.1f degs', subjectID, 0, 0));
    
    % Generate the OI from the wavefront map
    theOI = wvf2oiSpecial(theWVF, micronsPerDegree, pupilDiamMM);
    
    psfSupportWavelength = wavelengthsListToCompute;
end


function zCoeffs = zCoeffsForSubject(subjectID, whichEye, measPupilDiameterMM, subtractCentralRefraction, refractiveErrorMicrons)

    % Load the Z-coeffs for all subjects
    [~, ~, subject_coeffs] = wvfLoadThibosVirtualEyes(measPupilDiameterMM);
    
    if (strcmp(whichEye, ThibosOptics.constants.rightEye))
        zCoeffs = subject_coeffs.rightEye;
    else
        zCoeffs = subject_coeffs.leftEye;
    end
    
    zCoeffs = zCoeffs(:,subjectID);
    
    theDefocusZcoeffIndex = wvfOSAIndexToVectorIndex('defocus');

    if (subtractCentralRefraction)
        zCoeffs(theDefocusZcoeffIndex) = 0;
    end

    % Add refractive error
    if (abs(refractiveErrorMicrons)>0)
        if (subtractCentralRefraction)
            fprintf('Original defocus (central-refraction subtracted) (microns): %f, added defocus (microns): %f\n', zCoeffs(theDefocusZcoeffIndex), refractiveErrorMicrons);
        else
            fprintf('Original defocus (microns): %f, added defocus (microns): %f\n', zCoeffs(theDefocusZcoeffIndex), refractiveErrorMicrons);
        end
    end
    zCoeffs(theDefocusZcoeffIndex) = zCoeffs(theDefocusZcoeffIndex) + refractiveErrorMicrons;

   if (subjectID == 0)
        % Special case: diffraction-limited optics, possibly with some
        % refractive error
        % All zero coefficients
        zCoeffs = 0*zCoeffs;
        % Add refractive error
        zCoeffs(theDefocusZcoeffIndex) = refractiveErrorMicrons;
   end

  zCoeffs = zCoeffs(1:15);

end