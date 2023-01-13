function [PSFs, OTFs, xSfCyclesDeg, ySfCyclesDeg, xMinutes, yMinutes, theWVF] = ...
        computePSFandOTF(zcoeffsMicrons, wavelengthsListToCompute, wavefrontSpatialSamples, measPupilDiamMM, ...
        targetPupilDiamMM, measWavelength, showTranslation, varargin)
    
    p = inputParser;
    p.addParameter('doNotZeroCenterPSF', false, @islogical);
    p.addParameter('micronsPerDegree', 300, @isscalar);
    p.addParameter('flipPSFUpsideDown', false, @islogical);
    p.addParameter('rotatePSF90degs', false, @islogical);
    p.addParameter('upsampleFactor', [], @(x)(isempty(x) || ((isnumeric(x))&&(numel(x)==1)&&(x>0))));
    p.addParameter('noLCA',false,@islogical);
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
            'noLCA', noLCA);
    
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
        
        % Also set the PSF stored within theWVF
        theWVF.psf{wIndex} = theWavePSF;
    end
    
   
    xMinutes = xGridMinutes(1,:);
    yMinutes = yGridMinutes(:,1);
end