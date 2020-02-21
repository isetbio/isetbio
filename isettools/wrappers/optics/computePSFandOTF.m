function [PSFs, OTFs, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = computePSFandOTF(Zcoeffs, wavelengthsListToCompute, wavefrontSpatialSamples, measPupilDiamMM, targetPupilDiamMM, centeringWavelength, showTranslation)
    %% Compute WVF
    umPerDegree = 300;
    theWVF = makeWVF(wavefrontSpatialSamples, Zcoeffs, wavelengthsListToCompute, ...
            measPupilDiamMM, targetPupilDiamMM, umPerDegree, '');
    
    xSfCyclesPerRetinalMicron = wvfGet(theWVF, 'otf support', 'um', wavelengthsListToCompute(1));
    xSfCyclesDeg = xSfCyclesPerRetinalMicron * wvfGet(theWVF,'um per degree');
    ySfCyclesDeg = xSfCyclesDeg;
    [xSfGridCyclesDegGrid,ySfGridCyclesDegGrid] = meshgrid(xSfCyclesDeg, ySfCyclesDeg);
    
    if (~isempty(centeringWavelength))
        % Retrieve OTF at the centeringWavelength
        theCenteringOTF = wvfGet(theWVF, 'otf', centeringWavelength);
        theCenteringPSF = wvfGet(theWVF, 'psf', centeringWavelength);
        translationVector = []; showTranslation = false;
        [~, translationVector, ~, ~, ~] = otfWithZeroCenteredPSF(...
                    theCenteringOTF, theCenteringPSF, ...
                    translationVector, xSfGridCyclesDegGrid,ySfGridCyclesDegGrid, ...
                    showTranslation);
    else
        translationVector = [0 0];
    end
    
    for wIndex = 1:numel(wavelengthsListToCompute)
        theWaveOTF = wvfGet(theWVF, 'otf', wavelengthsListToCompute(wIndex));
        theWavePSF = wvfGet(theWVF, 'psf', wavelengthsListToCompute(wIndex));
        [theWaveOTF, ~, ~, ~,~] = ...
            otfWithZeroCenteredPSF(theWaveOTF, theWavePSF, translationVector,  xSfGridCyclesDegGrid, ySfGridCyclesDegGrid, showTranslation);
        
        [xGridMinutes,yGridMinutes, theWavePSF] = ...
            OtfToPsf(xSfGridCyclesDegGrid,ySfGridCyclesDegGrid,fftshift(theWaveOTF));
        if (wIndex == 1)
            OTFs = zeros(size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthsListToCompute));
            PSFs = zeros(size(theWavePSF,1), size(theWavePSF,2), numel(wavelengthsListToCompute));
        end
        
        OTFs(:,:,wIndex) = fftshift(theWaveOTF);
        PSFs(:,:,wIndex) = theWavePSF;
    end
    
    xMinutes = xGridMinutes(1,:);
    yMinutes = yGridMinutes(:,1);
end