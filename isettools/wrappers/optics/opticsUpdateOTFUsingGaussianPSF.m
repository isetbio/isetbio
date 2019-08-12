function optics = opticsUpdateOTFUsingGaussianPSF(optics, psfSigmaMicrons, maxSF, deltaSF, wavelengthSupport)
% Update the OTF of the passed optics struct with an OTF derived from a Gaussian PSF. 
%
% Syntax:
%   optics = opticsWithGaussianPSF(optics, psfSigmaMicrons, maxSF, deltaSF, wavelengthSupport)
%


    % Generate spatial frequency support
    N = maxSF/deltaSF;
    fList = unitFrequencyList(N);  % Should we add fList(end+1) = -fList(1) ?? (NPC)
    fList = fList * maxSF;
    [xSfGridCyclesDeg,ySfGridCyclesDeg] = meshgrid(fList, fList);
    fSupportCyclesPerDeg(:, :, 1) = xSfGridCyclesDeg;
    fSupportCyclesPerDeg(:, :, 2) = ySfGridCyclesDeg;
    
    % Compute the sigma of the PSF in minutes of arc, 1 deg = 60 minutes or arc
    psfSigmaMinutes = psfSigmaMicrons / optics.micronsPerDegree * 60;
    
    % Generate the 2D OTF (single wavelength) based on a Gaussian PSF
    theOTF = otfFromGaussianPSF(psfSigmaMinutes, fSupportCyclesPerDeg);
    
    % Reshape OTF in isetbio format
    theOTFIsetbioFormat = ifftshift(theOTF);
    otfData = zeros(size(theOTFIsetbioFormat,1), size(theOTFIsetbioFormat,2), numel(wavelengthSupport));
    for iW = 1:numel(wavelengthSupport)
        otfData(:,:,iW) = theOTFIsetbioFormat;
    end
    
    % Compute OTF frequency support in cycles/mm
    mmPerDegree = optics.micronsPerDegree * 1e-3;
    fSupportCyclesPerMM = fSupportCyclesPerDeg / mmPerDegree;
    fxCyclesPerMM = fSupportCyclesPerMM(1, :, 1);
    fyCyclesPerMM = fSupportCyclesPerMM(:, 1, 2);
    
    optics = opticsSet(optics, 'otf data',otfData);
    optics = opticsSet(optics, 'otffx', fxCyclesPerMM(:)');
    optics = opticsSet(optics, 'otffy', fyCyclesPerMM(:)');
    optics = opticsSet(optics, 'otfWave', wavelengthSupport);
end

function theOTF = otfFromGaussianPSF(psfSigmaMinutes, fSupportCyclesPerDeg)
    xSfGridCyclesDeg = fSupportCyclesPerDeg(:, :, 1);
    ySfGridCyclesDeg = fSupportCyclesPerDeg(:, :, 2);
    [xGridMinutes,yGridMinutes] = SfGridCyclesDegToPositionGridMinutes(xSfGridCyclesDeg,ySfGridCyclesDeg);
    thePSF = gaussianPSF(xGridMinutes, yGridMinutes, psfSigmaMinutes);
    [~,~,theOTF] = PsfToOtf(xGridMinutes,yGridMinutes,thePSF);
end

function thePSF = gaussianPSF(xGridMinutes, yGridMinutes, sigmaMinutes)
    radius = sqrt(xGridMinutes.^2 + yGridMinutes.^2);
    thePSF = normpdf(radius,0,sigmaMinutes);
    % Unit volume
    thePSF = thePSF/sum(thePSF(:));
end