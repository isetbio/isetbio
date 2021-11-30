function oi = oiFromPSF(psfEnsemble, wavelengthSupport, spatialSupportArcMin, ...
    pupilDiameterMM, umPerDegree)
% Create an oi from an ensemble of PSFs specified at a range of wavelenths
%
% Description:
%   Create an oi from an ensemble of PSFs specified at a range of wavelenths
%
% Inputs:
%   psfEnsemble             [nSamples x nSamples x nWaves] matrix of input PSFs
%   wavelengthSupport       vector of wavelengths at which the PSFs are measured 
%   spatialSupportArcMin    vector of spatial samples (in arc min) at which the PSFs are measured           
%   pupilDiameterMM         pupil diameter at which the PSFs are measured
%   umPerDegree             retinal magnification factor in microns/degree (~300 for the human, ~220 for the macaque)
%
% Outputs:
%  oi                       the generated oi

% History:
%   11/23/2021  Nicolas P. Cottaris, ISETBio Team. Wrote it.

    % Ensure that dimensions of input vectors are consistent
    assert(numel(spatialSupportArcMin) == size(psfEnsemble,1), ...
        sprintf('length of ''spatialSupportArcMin'' (%d) does not match rows of ''psfEnsemble'' (%d).', ...
                numel(spatialSupportArcMin), size(psfEnsemble,1)));
    assert(numel(spatialSupportArcMin) == size(psfEnsemble,2), ...
        sprintf('length of ''spatialSupportArcMin'' (%d) does not match cols of ''psfEnsemble'' (%d).', ...
                numel(spatialSupportArcMin), size(psfEnsemble,2)));
    assert(numel(wavelengthSupport) == size(psfEnsemble,3), ...
        sprintf('length of ''wavelengthSupport'' (%d) does not match z-dimension of ''psfEnsemble'' (%d).', ...
                numel(wavelengthSupport), size(psfEnsemble,3)));
            
    % Generate position grids
    [xGridMinutes,yGridMinutes] = meshgrid(spatialSupportArcMin,spatialSupportArcMin);
    
    for wIndex = 1:numel(wavelengthSupport)
        % Retrieve PSF at this wavelength
        theWavePSF = psfEnsemble(:,:,wIndex);
        
        % Convert PSF to OTF
        [xSfGridCyclesDeg, ySfGridCyclesDeg, theWaveOTF] = ...
            PsfToOtf(xGridMinutes, yGridMinutes, theWavePSF/sum(theWavePSF(:)));
        
        % Allocate memory
        if (wIndex == 1)
            otf = zeros(size(theWaveOTF,1), size(theWaveOTF,2), numel(wavelengthSupport));
        end
        
        % Isetbio wants the otf with (0, 0) sf at the upper left.  We
        % accomplish this by applying ifftshift to the wvf centered format.
        otf(:, :, wIndex) = ifftshift(theWaveOTF);    
    end
    
    % Compute spatial frequency support in cycles/mm
    xSfCyclesDeg = xSfGridCyclesDeg(1,:);
    ySfCyclesDeg = ySfGridCyclesDeg(:,1);
    spatialFrequencySupportXcyclesPerMM = xSfCyclesDeg / (1e-3*umPerDegree);
    spatialFrequencySupportYcyclesPerMM = ySfCyclesDeg / (1e-3*umPerDegree);
    
    % Create oi
    oi = oiCreate;
    oi = oiSet(oi, 'name', 'oi from PSF');

    % Set the OTF data and support parameters.
    oi = oiSet(oi, 'optics OTF fx', spatialFrequencySupportXcyclesPerMM);
    oi = oiSet(oi, 'optics OTF fy', spatialFrequencySupportYcyclesPerMM);
    oi = oiSet(oi, 'optics OTF wave', wavelengthSupport);
    oi = oiSet(oi, 'optics otfdata', otf);
    
    % Set the oi wavelength support
    oi = oiSet(oi, 'wave', wavelengthSupport);
    
    % Set the oi fNumber 
    focalLengthMM = (umPerDegree * 1e-3) / (2 * tand(0.5));
    fNumber = focalLengthMM/pupilDiameterMM;
    oi = oiSet(oi, 'optics fnumber', fNumber);
    
    % Set the optics fNumber
    optics = oiGet(oi, 'optics');
    optics = opticsSet(optics, 'fnumber', fNumber);
    
    % Set the optics focal length (in meters)
    optics = opticsSet(optics, 'focalLength', focalLengthMM*1e-3);
    oi = oiSet(oi, 'optics', optics);
end


