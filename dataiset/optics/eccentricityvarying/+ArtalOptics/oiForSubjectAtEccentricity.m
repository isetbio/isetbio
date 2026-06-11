function [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength, zCoeffs] = oiForSubjectAtEccentricity(subjectID, whichEye, ecc, ...
    pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin)

    % Parse input
    p = inputParser;
    p.addRequired('subjectID', @(x)(isscalar(x)&&(x>=0)&&(x<=130)));
    p.addRequired('whichEye', @(x)(ischar(x)&&(ismember(x,{ArtalOptics.constants.rightEye, ArtalOptics.constants.leftEye}))));
    p.addRequired('ecc', @(x)(isnumeric(x)&&((numel(x) == 1)||(numel(x) == 2))));
    p.addRequired('pupilDiamMM', @(x)(isscalar(x)&&(x>=1)&&(x<=4)));
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