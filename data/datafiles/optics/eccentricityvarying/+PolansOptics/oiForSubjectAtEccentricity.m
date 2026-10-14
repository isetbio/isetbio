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
    p.addParameter('subtractCentralRefraction3D', false, @islogical);
    p.addParameter('withZeroedPistonAndTiltZernikeCoefficients', false, @islogical);
    p.addParameter('zeroCenterPSF', false, @islogical);
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
    p.addParameter('noLCA', false, @islogical);
    p.addParameter('refractiveErrorDiopters', 0, @isnumeric);
    p.addParameter('obliqueAndVerticalAstigmatismErrorsDiopters', [0 0], @isnumeric);
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
    subtractCentralRefraction3D = p.Results.subtractCentralRefraction3D;
    
    zeroCenterPSF = p.Results.zeroCenterPSF;
    withZeroedPistonAndTiltZernikeCoefficients = p.Results.withZeroedPistonAndTiltZernikeCoefficients;
    flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
    upsampleFactor = p.Results.upsampleFactor;
    noLCA = p.Results.noLCA;

    % Refractive (defocus) error, if desired
    defocusErrorDiopters = p.Results.refractiveErrorDiopters;
    defocusErrorMicrons = wvfDefocusDioptersToMicrons(defocusErrorDiopters, pupilDiamMM);
    
    % Astigmatic errors, if desired
    obliqueAndVerticalAstigmatismErrorsDiopters = p.Results.obliqueAndVerticalAstigmatismErrorsDiopters;
    obliqueAndVerticalAstigmatismErrorsMicrons = wvfAstigmatismDioptersToMicrons(obliqueAndVerticalAstigmatismErrorsDiopters, pupilDiamMM);

    % Obtain z-coeffs at desired eccentricity
    zCoeffs = zCoeffsForSubjectAtEcc(subjectID, ecc, ...
        subtractCentralRefraction, subtractCentralRefraction3D, ...
        defocusErrorMicrons, obliqueAndVerticalAstigmatismErrorsMicrons, pupilDiamMM);
    if (subjectID == 0)
        measurementPupilDiameterMM = pupilDiamMM;
    else
        measurementPupilDiameterMM = PolansOptics.constants.measurementPupilDiamMM;
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
             'withZeroedPistonAndTiltZernikeCoefficients', withZeroedPistonAndTiltZernikeCoefficients, ...
             'micronsPerDegree', micronsPerDegree, ...
             'flipPSFUpsideDown', flipPSFUpsideDown, ...
             'upsampleFactor', upsampleFactor, ...
             'noLCA', noLCA, ...
             'name', sprintf('Polans subject %d, eccentricity: %2.1f,%2.1f degs', subjectID, ecc(1), ecc(2)));

    % Need to put defocus back in for subject 0 case, because it got 
    % zeroed out above.  Easiest to do that 
    if (subjectID == 0)
        theWVF = wvfSet(theWVF,'zcoeffs', defocusErrorMicrons, {'defocus'});
    end
    
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


function  interpolatedZcoeffs = zCoeffsForSubjectAtEcc(subjectID, ecc, ...
    subtractCentralRefraction, subtractCentralRefraction3D, ...
    defocusErrorMicrons, obliqueAndVerticalAstigmatismErrorsMicrons, pupilDiamMM)

    % Get original z-coeffs at all measured eccentricities
    if (subjectID == 0)
        [zMap, zCoeffIndices] = PolansOptics.constants.ZernikeCoefficientsMap(1);
    else
        [zMap, zCoeffIndices] = PolansOptics.constants.ZernikeCoefficientsMap(subjectID);
    end
    
    % Check that the coefficients do not contain NaNs
    atLeastOneNaNCoefficient = any(isnan(zMap(:)));
    
    % Check that the coefficients are not all 0
    allZeroCoefficients = (sum(abs(zMap(:))) == 0);
    
    if (atLeastOneNaNCoefficient || allZeroCoefficients)
        interpolatedZcoeffs = [];
        return;
    end

    zCoeffsNum = size(zMap,3);

    % Interpolate zMap at desired ecc
    [X,Y] = meshgrid(...
        PolansOptics.constants.measurementHorizontalEccentricities, ...
        PolansOptics.constants.measurementVerticalEccentricities);

    % index of coeffs at ecc = (0,0)
    indexOfZeroEcc = (X==0) & (Y==0);

    interpolatedZcoeffs = zeros(1, 30);

    % The Polans Z-coeff indices start from 0, not 1
    theDefocusZcoeffIndex = wvfOSAIndexToVectorIndex('defocus') - 1;
    theObliqueAstigmatismZcoeffIndex = wvfOSAIndexToVectorIndex('oblique_astigmatism') - 1;
    theVerticalAstigmatismZcoeffIndex = wvfOSAIndexToVectorIndex('vertical_astigmatism') - 1;

    for zIndex = 1:zCoeffsNum
         % Retrieve the XY spatial map for this z-coeff
         z2Dmap = squeeze(zMap(:,:,zIndex));
         
         % Subtract central defocus from all spatial positions
         if ((zCoeffIndices(zIndex) == theDefocusZcoeffIndex) && (subtractCentralRefraction||subtractCentralRefraction3D))
             z2Dmap = z2Dmap - z2Dmap(indexOfZeroEcc);
         end
         
         % Subtract oblique astigmatism from all spatial positions
        if ((zCoeffIndices(zIndex) == theObliqueAstigmatismZcoeffIndex) && (subtractCentralRefraction3D))
            %fprintf('Subtracting central oblique astigmatism: %2.2f D (%d/%d)\n', wvfAstigmatismMicronsToDiopters(z1Dmap(indexOfZeroEcc), pupilDiamMM), theObliqueAstigmatismZcoeffIndex, zIndex);
            z2Dmap = z2Dmap - z2Dmap(indexOfZeroEcc);
        end

        % Subtract vertical astigmatism from all spatial positions
        if ((zCoeffIndices(zIndex) == theVerticalAstigmatismZcoeffIndex) && (subtractCentralRefraction3D))
            %fprintf('Subtracting central vertical astigmatism: %2.2f D (%d/%d)\n', wvfAstigmatismMicronsToDiopters(z1Dmap(indexOfZeroEcc), pupilDiamMM),  theVerticalAstigmatismZcoeffIndex, zIndex);
            z2Dmap = z2Dmap - z2Dmap(indexOfZeroEcc);
        end

         % Add defocus error
         if (zCoeffIndices(zIndex) == theDefocusZcoeffIndex)
             if (abs(defocusErrorMicrons)>0)
                %if (subtractCentralRefraction)
                %    fprintf('Original defocus at eccXY = (0,0) (central-refraction subtracted) (microns): %f, added defocus (microns): %f\n', z2Dmap(indexOfZeroEcc), defocusErrorMicrons);
                %else
                %    fprintf('Original defocus at eccXY = (0,0) (microns): %f, added defocus (microns): %f\n', z2Dmap(indexOfZeroEcc), defocusErrorMicrons);
                %end
             end
             z2Dmap = z2Dmap + defocusErrorMicrons;
         end
         

        % Add oblique astigmatism error
        if (zCoeffIndices(zIndex) == theObliqueAstigmatismZcoeffIndex) && (abs(obliqueAndVerticalAstigmatismErrorsMicrons(1))>0)
            %fprintf(2,'Original oblique astigmatism: %f (microns), added (microns): %f\n', z2Dmap(indexOfZeroEcc), obliqueAndVerticalAstigmatismErrorsMicrons(1));
            z2Dmap = z2Dmap + obliqueAndVerticalAstigmatismErrorsMicrons(1);
        end

        % Add vetical astigmatism error
        if (zCoeffIndices(zIndex) == theVerticalAstigmatismZcoeffIndex) && (abs(obliqueAndVerticalAstigmatismErrorsMicrons(2))>0)
            %fprintf(2,'Original vertical astigmatism: %f (microns), added (microns): %f\n', z2Dmap(indexOfZeroEcc), obliqueAndVerticalAstigmatismErrorsMicrons(2));
            z2Dmap = z2Dmap + obliqueAndVerticalAstigmatismErrorsMicrons(2);
        end

        % Interpolate the XY map at the desired eccentricity.
        % Note the + 1 added to the z-coeff. This is to ensure that
        % we are addressing the correct indices expected by the wvf object
        interpolatedZcoeffs(zCoeffIndices(zIndex)+1) = interp2(X,Y,z2Dmap, ecc(1), ecc(2));
    end
    
    if (subjectID == 0)
        % All zero coefficients
        interpolatedZcoeffs = 0*interpolatedZcoeffs;

        % Add defocus error
        interpolatedZcoeffs(wvfOSAIndexToVectorIndex('defocus')) = defocusErrorMicrons;

        % Add oblique_astigmatism  error
        interpolatedZcoeffs(wvfOSAIndexToVectorIndex('oblique_astigmatism')) = obliqueAndVerticalAstigmatismErrorsMicrons(1);

        % Add vertical_astigmatism  error
        interpolatedZcoeffs(wvfOSAIndexToVectorIndex('vertical_astigmatism')) = obliqueAndVerticalAstigmatismErrorsMicrons(2);

    end
end