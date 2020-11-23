function [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY] = oiForSubjectAtEccentricity(subjectID, ecc, ...
    pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin)

    % Parse input
    p = inputParser;
    p.addRequired('subjectID', @(x)(isscalar(x)&&(x>=1)&&(x<=10)));
    p.addRequired('ecc', @(x)(isnumeric(x)&&(numel(x) == 2)));
    p.addRequired('pupilDiamMM', @(x)(isscalar(x)&&(x>=1)&&(x<=4)));
    p.addRequired('wavelengthsListToCompute', @(x)(isnumeric(x)));
    p.addRequired('micronsPerDegree', @(x)(isscalar(x)));
    p.addParameter('inFocusWavelength', 550, @isscalar);
    p.addParameter('wavefrontSpatialSamples', 801, @isscalar)
    p.addParameter('subtractCentralRefraction', true, @islogical);
    p.addParameter('noLCA', false, @islogical);
    p.parse(subjectID, ecc, pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin{:});
    
    inFocusWavelength = p.Results.inFocusWavelength;
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
    subtractCentralRefraction = p.Results.subtractCentralRefraction;
    noLCA = p.Results.noLCA;
   
    
    % Obtain z-coeffs at desired eccentricity
    zCoeffs = zCoeffsForSubjectAtEcc(subjectID, ecc, subtractCentralRefraction);
    
    % Compute PSF and WVF from z-Coeffs for the desired pupil and wavelenghts
    [thePSF, ~, ~,~, psfSupportMinutesX, psfSupportMinutesY, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             PolansOptics.constants.measurementPupilDiamMM, ...
             pupilDiamMM, inFocusWavelength, false, ...
             'doNotZeroCenterPSF', true, ...
             'micronsPerDegree', micronsPerDegree);
    
    % Remove wavelength-dependent defocus if noLCA is set
    if (noLCA)
        % Set all PSFs to the PSF at the in-focus wavelenth
        [~,wTarget] = min(abs(wavelengthsListToCompute-inFocusWavelength));
        targetPSF = thePSF(:,:,wTarget);
        for waveIndex = 1:numel(wavelengthsListToCompute)
            theWVF.psf{waveIndex} = targetPSF;
            thePSF(:,:,waveIndex) = targetPSF;
        end
    end
    
    % Generate the OI from the wavefront map
    theOI = wvf2oi(theWVF);
end

function  interpolatedZcoeffs = zCoeffsForSubjectAtEcc(subjectID, ecc, subtractCentralRefraction)

    % Get original z-coeffs at all measured eccentricities
    [zMap, zCoeffIndices] = PolansOptics.constants.ZernikeCoefficientsMap(subjectID);
    zCoeffsNum = size(zMap,3);
    
    % Interpolate zMap at desired ecc
    [X,Y] = meshgrid(...
        PolansOptics.constants.measurementHorizontalEccentricities, ...
        PolansOptics.constants.measurementVerticalEccentricities);
    
    interpolatedZcoeffs = zeros(1, 30);
    for zIndex = 1:zCoeffsNum
         % Retrieve the XY map for this z-coeff
         zz = squeeze(zMap(:,:,zIndex));
         % The 4-th z-coeff is defocus. Subtract central defocus from all
         % spatial positions
         if ((zCoeffIndices(zIndex) == 4) && (subtractCentralRefraction))
             idx = find((X==0) & (Y==0));
             zz = zz - zz(idx);
         end
         % Interpolate the XY map at the desired eccentricity
         interpolatedZcoeffs(zCoeffIndices(zIndex)+1) = interp2(X,Y,zz, ecc(1), ecc(2));
     end
end