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
    p.addParameter('micronsPerDegree', 300, @isscalar);
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('rotatePSF90degs', false, @islogical);
    p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
    p.addParameter('noLCA',false,@islogical);  % This should become human lca
    p.addParameter('name', 'noname', @ischar);
    p.parse(varargin{:});
    doNotZeroCenterPSF = p.Results.doNotZeroCenterPSF;
    flipPSFUpsideDown = p.Results.flipPSFUpsideDown;
    rotatePSF90degs = p.Results.rotatePSF90degs;
    upsampleFactor = p.Results.upsampleFactor;
    umPerDegree = p.Results.micronsPerDegree;
    noLCA = p.Results.noLCA;
    name = p.Results.name;
    
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

            % In the ISETCam branch, wvfGet returns the OTF with zero sf in the
            % middle.  But this routine was expecting it to have zero in the
            % upper left, as accomplished by ifftshift. Thus we apply an
            % ifftshift() here.
            theCenteringOTF = ifftshift(wvfGet(theWVF, 'otf', theWavelengthSupport(idx)));
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
        % In the ISETCam branch, wvfGet returns the OTF with zero sf in the
        % middle.  But this routine was expecting it to have zero in the
        % upper left, as accomplished by ifftshift. Thus we apply an
        % ifftshift here.
        theWaveOTF = ifftshift(wvfGet(theWVF, 'otf', wavelengthsListToCompute(wIndex)));

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
        
        % The fftshift() puts the OTF back with zero sf in the middle.
        % I'm not sure why this was being done in the ISETBio version, but
        % it seems right. And, OTF was never picked up on return by any of
        % the existing calls to this routine, so I think it is OK to have
        % it now match the format returned by wvfGet().
        OTFs(:,:,wIndex) = fftshift(theWaveOTF);
        PSFs(:,:,wIndex) = theWavePSF;
        
        % Also set the PSF stored within theWVF
        theWVF.psf{wIndex} = theWavePSF;
    end
    
   
    xMinutes = xGridMinutes(1,:);
    yMinutes = yGridMinutes(:,1);
end