function [PSFs, OTFs, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = ...
        computePSFandOTF(zcoeffsMicrons, wavelengthsListToCompute, wavefrontSpatialSamples, measPupilDiamMM, ...
        targetPupilDiamMM, measWavelength, showTranslation, varargin)
% This function centers the PSF.  It appears to do so by adjusting the
% phase of the OTF, so the OTF an PSF remain matched.
%
% When we changed from ISETBio to ISETCam, the convention of the OTF
% returned by wvfGet changed from zero sf in the upper left (ISETBio) to
% zero sf in the middle (ISETCam).  That broke this routine, but we have
% fixed it as commented below by putting a call to ifftshift() around the 
% OTFs returned by wvfGet.  We are currently thinking this is all that is
% needed.  Comments about those places are below.

    p = inputParser;
    p.addParameter('doNotZeroCenterPSF', false, @islogical);
    p.addParameter('withZeroedPistonAndTiltZernikeCoefficients', false, @islogical);
    p.addParameter('micronsPerDegree', 300, @isscalar);
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('rotatePSF90degs', false, @islogical);
    p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
    p.addParameter('noLCA',false,@islogical);  % This should become human lca
    p.addParameter('name', 'noname', @ischar);
    p.parse(varargin{:});
    withZeroedPistonAndTiltZernikeCoefficients = p.Results.withZeroedPistonAndTiltZernikeCoefficients;
    doNotZeroCenterPSF = p.Results.doNotZeroCenterPSF;
    flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
    rotatePSF90degs = p.Results.rotatePSF90degs;
    upsampleFactor = p.Results.upsampleFactor;
    umPerDegree = p.Results.micronsPerDegree;
    noLCA = p.Results.noLCA;
    name = p.Results.name;
    

    % Zero the first 3 Zernike coeffs, and set the  doNotZeroCenterPSF to false
    if (withZeroedPistonAndTiltZernikeCoefficients)
        % Zero the first 3 Zernike coefficients: 
        % piston (1), 
        % vertical tilt (2), and 
        % horizontal tilt (3)
        zcoeffsMicrons(1:3) = 0;
        if (doNotZeroCenterPSF)
            % User has passed ''doNotZeroCenterPSF'' = true, and
            % ''withZeroedPistonAndTiltZernikeCoefficients'' = true
            % Notify the user that we are overriding the value of '''doNotZeroCenterPSF''
            fprintf('\n\ncomputePSFandOTF(): overiding the passed value ''doNotZeroCenterPSF'' from true to false, because the passed ''withZeroedPistonAndTiltZernikeCoefficients'' is set to true.\n\n');
            doNotZeroCenterPSF = false;
        end
    else
        if (doNotZeroCenterPSF == false)
            % User wants zero-centered PSF. Suggest switching to setting 
            % ''withZeroedPistonAndTiltZernikeCoefficients'' to true
            % instead.
            fprintf('\n\ncomputePSFandOTF(): the passed value ''doNotZeroCenterPSF'' is set to false and ''withZeroedPistonAndTiltZernikeCoefficients'' is also false.\n Consider setting ''withZeroedPistonAndTiltZernikeCoefficients'' to true instead.\n\n');
        else
            % Old behavior. Notidy user that
            % ''withZeroedPistonAndTiltZernikeCoefficients''  will default to true in the future.
            fprintf('\n\ncomputePSFandOTF(): ''withZeroedPistonAndTiltZernikeCoefficients'' still defaults to false, but soon will default to true.\n\n');
        end
    end

    %% Compute WVF
    theWVF = makeWVF(wavefrontSpatialSamples, zcoeffsMicrons, measWavelength, wavelengthsListToCompute, ...
            measPupilDiamMM, targetPupilDiamMM, umPerDegree, name, ...
            'upsampleFactor', upsampleFactor, ...
            'flipPSFUpsideDown', flipPSFUpsideDown, ...
            'rotatePSF90degs', rotatePSF90degs , ...
            'humanlca', ~noLCA);
    
    xSfCyclesPerRetinalMicron = wvfGet(theWVF, 'otf support', 'um', wavelengthsListToCompute(1));
    xSfCyclesDeg = xSfCyclesPerRetinalMicron * wvfGet(theWVF,'um per degree');
    ySfCyclesDeg = xSfCyclesDeg;
    [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfCyclesDeg, ySfCyclesDeg);
    

    if (doNotZeroCenterPSF)
        translationVector = [0 0];
    else
        if (~isempty(measWavelength))
            % Retrieve OTF closest to the measurement wavelength]
            theWavelengthSupport = wvfGet(theWVF, 'wave');
            [~,idx] = min(abs(theWavelengthSupport-measWavelength));

            % 04/02/24, DHB: Remove ifftshift to match change to wvfGet (where
            % the ifftshift went back in)
            theCenteringOTF = wvfGet(theWVF, 'otf', theWavelengthSupport(idx));
            theCenteringPSF = wvfGet(theWVF, 'psf', theWavelengthSupport(idx));
            translationVector = []; showTranslation = false;
            [~, translationVector, ~, ~, ~] = otfWithZeroCenteredPSF(...
                theCenteringOTF, theCenteringPSF, ...
                translationVector, xSfGridCyclesDegGrid,ySfGridCyclesDegGrid, ...
                showTranslation);
        else
            translationVector = [0 0];
        end
    end
    
    for wIndex = 1:numel(wavelengthsListToCompute)
        % 04/02/24, DHB: Remove ifftshift to match change to wvfGet (where
        % the ifftshift went back in)
        theWaveOTF = wvfGet(theWVF, 'otf', wavelengthsListToCompute(wIndex));

        % The rest of this code then works as it did in the ISETBio branch,
        % where the OTF had zero sf in the upper left.
        theWavePSF = wvfGet(theWVF, 'psf', wavelengthsListToCompute(wIndex));
        [theWaveOTF, ~, ~, ~,~] = ...
            otfWithZeroCenteredPSF(theWaveOTF, theWavePSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation);
        
        % OtfToPsf expects zero sf in the middle, hence the fftshift() in
        % this call, as was there in the ISETBio version.
        [xGridMinutes,yGridMinutes, theWavePSF] = ...
            OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));
        if (wIndex == 1)
            OTFs = zeros(size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthsListToCompute));
            PSFs = zeros(size(theWavePSF,1), size(theWavePSF,2), numel(wavelengthsListToCompute));
        end
        
        % 04/02/24: Take out fftshift on OTF, because we now are putting DC
        % at 1,1
        OTFs(:,:,wIndex) = theWaveOTF;
        PSFs(:,:,wIndex) = theWavePSF;
        
        % Also set the PSF stored within theWVF
        theWVF.psf{wIndex} = theWavePSF;
    end
    
   
    xMinutes = xGridMinutes(1,:);
    yMinutes = yGridMinutes(:,1);
end