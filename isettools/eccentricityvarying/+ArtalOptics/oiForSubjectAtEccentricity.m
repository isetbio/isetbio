function [theOI, thePSF, psfSupportMinutesX, psfSupportMinutesY, psfSupportWavelength] = oiForSubjectAtEccentricity(subjectID, whichEye, ecc, ...
    pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin)

    % Parse input
    p = inputParser;
    p.addRequired('subjectID', @(x)(isscalar(x)&&(x>=1)&&(x<=130)));
    p.addRequired('whichEye', @(x)(ischar(x)&&(ismember(x,{ArtalOptics.constants.rightEye, ArtalOptics.constants.leftEye}))));
    p.addRequired('ecc', @(x)(isnumeric(x)&&(numel(x) == 2)));
    p.addRequired('pupilDiamMM', @(x)(isscalar(x)&&(x>=1)&&(x<=4)));
    p.addRequired('wavelengthsListToCompute', @(x)(isnumeric(x)));
    p.addRequired('micronsPerDegree', @(x)(isscalar(x)));
    p.addParameter('inFocusWavelength', 550, @isscalar);
    p.addParameter('wavefrontSpatialSamples', 801, @isscalar)
    p.addParameter('subtractCentralRefraction', true, @islogical);
    p.addParameter('zeroCenterPSF', false, @islogical);
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('noLCA', false, @islogical);
    p.parse(subjectID, whichEye, ecc, pupilDiamMM, wavelengthsListToCompute, micronsPerDegree, varargin{:});
    
    inFocusWavelength = p.Results.inFocusWavelength;
    wavefrontSpatialSamples = p.Results.wavefrontSpatialSamples;
    subtractCentralRefraction = p.Results.subtractCentralRefraction;
    zeroCenterPSF = p.Results.zeroCenterPSF;
    flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
    noLCA = p.Results.noLCA;
    
    
    % Flip vertical eccentricity
    % Obtain z-coeffs at desired eccentricity
    zCoeffs = zCoeffsForSubjectAtEcc(subjectID, whichEye, ecc(1), subtractCentralRefraction);
    
    % Compute PSF and WVF from z-Coeffs for the desired pupil and wavelenghts
    [thePSF, ~, ~,~, psfSupportMinutesX, psfSupportMinutesY, theWVF] = ...
        computePSFandOTF(zCoeffs, ...
             wavelengthsListToCompute, wavefrontSpatialSamples, ...
             ArtalOptics.constants.measurementPupilDiamMM, ...
             pupilDiamMM, inFocusWavelength, false, ...
             'doNotZeroCenterPSF', ~zeroCenterPSF, ...
             'micronsPerDegree', micronsPerDegree, ...
             'flipPSFUpsideDown', flipPSFUpsideDown, ...
             'name', sprintf('Artal subject %d, eccentricity: %2.1f,%2.1f degs', subjectID, ecc(1), ecc(2)));
    
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
    theOI = wvf2oiSpecial(theWVF, micronsPerDegree, pupilDiamMM);
    
    psfSupportWavelength = wavelengthsListToCompute;
end

function theOI = wvf2oiSpecial(theWVF, umPerDegree, pupilDiameterMM)

    % Generate oi from the wvf
    theOI = wvf2oi(theWVF);
    
    % Adjust the OI's fNumber and focalLength to be consistent with the
    % micronsPerDegree and pupilDiameter of the WVF
    optics = oiGet(theOI, 'optics');
    focalLengthMM = (umPerDegree * 1e-3) / (2 * tand(0.5));
    fLengthMeters = focalLengthMM * 1e-3;
    pupilRadiusMeters = (pupilDiameterMM / 2) * 1e-3;
    pupilDiameterMeters = 2 * pupilRadiusMeters;
    optics = opticsSet(optics, 'fnumber', fLengthMeters / pupilDiameterMeters);
    optics = opticsSet(optics, 'focalLength', fLengthMeters);
    theOI = oiSet(theOI, 'optics', optics);

    %heightDegs = oiGet(theOI, 'hangular');
    %heightMicrons = oiGet(theOI, 'height')*1e6;
    %fprintf('Achievend microns per deg: %f (desired: %f)\n', heightMicrons/heightDegs, umPerDegree);  
end


function  interpolatedZcoeffs = zCoeffsForSubjectAtEcc(subjectID, whichEye, ecc, subtractCentralRefraction)

    % Get original z-coeffs at all measured eccentricities
    [zMap, zCoeffIndices] = ArtalOptics.constants.ZernikeCoefficientsMap(subjectID, whichEye);
    zCoeffsNum = size(zMap,2);
    
    % Index of coeffs at ecc = (0,0)
    indexOfZeroEcc = find(ArtalOptics.constants.measurementHorizontalEccentricities==0);
    
    interpolatedZcoeffs = zeros(1, 30);
    for zIndex = 1:zCoeffsNum
         % Retrieve the XY spatial map for this z-coeff
         z2Dmap = squeeze(zMap(:,zIndex));
         
         % The 4-th z-coeff is defocus. Subtract central defocus from all
         % spatial positions
         if ((zCoeffIndices(zIndex) == 4) && (subtractCentralRefraction))
             z2Dmap = z2Dmap - z2Dmap(indexOfZeroEcc);
         end
         % Interpolate the XY map at the desired eccentricity.
         interpolatedZcoeffs(zCoeffIndices(zIndex)+1) = interp1(ArtalOptics.constants.measurementHorizontalEccentricities,z2Dmap, ecc(1));
    end
     
end