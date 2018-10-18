function t_diffractionAcrossWavelengthPsfAndOtf
% Test different ways of getting diffraction limited optics with defocus
%
% Syntax:
%   t_diffractionAcrossWavelengthPsfAndOtf
%
% Description:
%    Shows how to get diffraction limited optics plus the effect of defocus
%    across wavelengths.
%
%    Also illustrates that we can get the same answer after inserting the
%    underlying wavefront optics structure into an oi and getting it back
%    out again, as well as go back and forth between psf and otf.
%
%    This was originally a test and helped us shake out a few numerical
%    bugs in our conversions.
%
%    The first two rows (otfs) of the big plot figure should match, as
%    should the next four rows (psfs).
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

% History:
%    01/25/18  ncp  Wrote it.
%    01/26/18  dhb  Finish up.
%              dhb  Remove dependence on NicePlot so we can include in
%                   isetbio.
%    09/13/18  jnm  Formatting

    % Close 
    close all
    
    % Set up figure
    hFig = figure(1);
    clf;
    set(hFig, 'Position', [10 10 1400 1600]);
    
    % Wavelengths and other parameters
    targetWavelengths = [450 470 550 590  610];
    spatialSamples = 501;
    
    for wIndex = 1:numel(targetWavelengths)
        test(targetWavelengths(wIndex), spatialSamples, wIndex, ...
            numel(targetWavelengths))
    end
end
 
function test(targetWavelength, spatialSamples, col, cols)
% Create wvf struct for diffraction limited optics at specified wavelength
%
% Syntax:
%   test(targetWavelength, spatialSamples, col, cols)
%
% Description:
%    A support function for t_diffractionAcrossWavelengthPsfAndOtf that
%    will create a WVF structure for diffraction limited optics at the
%    specified wavelength.
%
% Inputs:
%    targetWavelength - Numeric. The target wavelength, in nanometers.
%    spatialSamples   - Numeric. The number of spatial samples.
%    col              - Numeric. Which column of the figure to display in.
%    cols             - Numeric. The total number of columns in display.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%

    theWVF = wvfCreate('calc wavelengths', targetWavelength, ...
        'spatial samples', spatialSamples);
    theWVF = wvfComputePSF(theWVF);
    
    % Get the PSF and the spatial sampling for it from the wvf structure.
    % Note that spatial sampling can in principle be wavelength dependent.
    psfTargetWFromWvfPsf = wvfGet(theWVF, 'psf', targetWavelength);
    psfMinutesFromWvfPsf = ...
        wvfGet(theWVF, 'psf angular samples', 'min', targetWavelength);
    [xGridMinutesFromWvfPsf, yGridMinutesFromWvfPsf] = ...
        meshgrid(psfMinutesFromWvfPsf, psfMinutesFromWvfPsf);
    
    % Get OTF from the WVF together with the OTF sampling at the target
    % wavelength. Note, otf sampling can in principle be wavelength
    % dependent.
    otfTargetWFromWvf = wvfGet(theWVF, 'otf', targetWavelength);
    xSfCyclesDegFromWvfOtf = ...
        wvfGet(theWVF, 'otf support', 'um', targetWavelength) ...
        * wvfGet(theWVF, 'um per degree');
    ySfCyclesDegFromWvfOtf = xSfCyclesDegFromWvfOtf;
    [xSfGridCyclesDegFromWvfOtf, ySfGridCyclesDegFromWvfOtf] = ...
        meshgrid(xSfCyclesDegFromWvfOtf, ySfCyclesDegFromWvfOtf);
    [xGridMinutesFromWvfOtf, yGridMinutesFromWvfOtf, ...
        psfTargetWFromWvfOtf] = OtfToPsf(xSfGridCyclesDegFromWvfOtf, ...
        ySfGridCyclesDegFromWvfOtf, fftshift(otfTargetWFromWvf));
    
    % Check. The otf we get from the wvf should be the same thing we get by
    % applying PsfToOtf to the wvf psf, followed by ifftshift.
    [~, ~, otfTargetWFromWvfCheck] = ...
        PsfToOtf([], [], psfTargetWFromWvfPsf);
    otfTargetWFromWvfCheck = ifftshift(otfTargetWFromWvfCheck);
    if (any(otfTargetWFromWvfCheck(:) ~= otfTargetWFromWvf(:)))
        error('Mystery difference in otfs\n');
    end
    
    % And if the two OTFs match, then we should be able to go from the otf
    % back to the PSF we started with.
    [~, ~, psfTargetWFromWvfPsfCheck] = ...
        OtfToPsf([], [], fftshift(otfTargetWFromWvf));
    if (max(abs(psfTargetWFromWvfPsfCheck(:) - ...
            psfTargetWFromWvfPsf(:))) > 1e-10)
        error('Failure of OtfToPsf/PsfToOtf to self-invert\n');
    end
    
    % Convert the otf we got from the wvf to the psf directly.
    [xSfGridCyclesDegFromWvfOtf, ySfGridCyclesDegFromWvfOtf] = ...
        meshgrid(xSfCyclesDegFromWvfOtf, ySfCyclesDegFromWvfOtf);
    [xGridMinutesFromWvfOtf, yGridMinutesFromWvfOtf, ...
        psfTargetWFromWvfOtf] = OtfToPsf(xSfGridCyclesDegFromWvfOtf, ...
        ySfGridCyclesDegFromWvfOtf, fftshift(otfTargetWFromWvf));
    
    % Convert wvf to oi and then snag the optics structure out
    optics = oiGet(wvf2oi(theWVF), 'optics');
    
    % Get psf from optics structure
    psfTargetWFromOpticsPsf = ...
        opticsGet(optics, 'psf data', targetWavelength);
    psfDegreesFromOpticsPsf = opticsGet(optics, 'psf support', 'um');
    xGridMinutesFromOpticsPsf = 60 * psfDegreesFromOpticsPsf{1} ...
        / wvfGet(theWVF, 'um per degree');
    yGridMinutesFromOpticsPsf = 60 * psfDegreesFromOpticsPsf{2} ...
        / wvfGet(theWVF, 'um per degree');
    
    % Get OTF and support out of the optics structure
    otfTargetWFromOptics = opticsGet(optics, 'otf data', targetWavelength);
    xSfCyclesDegFromOpticsOtf = ...
        opticsGet(optics, 'otf fx', 'cyclesperdeg');
    ySfCyclesDegfromOpticsOtf = ...
        opticsGet(optics, 'otf fy', 'cyclesperdeg');
 
    % Convert the otf we got from the optics to the psf directly.
    [xSfGridCyclesDegFromOpticsOtf, ySfGridCyclesDegFromOpticsOtf] = ...
        meshgrid(xSfCyclesDegFromOpticsOtf, ySfCyclesDegfromOpticsOtf);
    [xGridMinutesFromOpticsOtf, yGridMinutesFromOpticsOtf, ...
        psfTargetWFromOpticsOtf] = OtfToPsf(...
        xSfGridCyclesDegFromOpticsOtf, ySfGridCyclesDegFromOpticsOtf, ...
        fftshift(otfTargetWFromOptics));
    
    % Common plot normalizer for psf plots
    maxAll = max([max(psfTargetWFromWvfPsf(:)), ...
        max(psfTargetWFromWvfOtf(:)), ...
        max(psfTargetWFromOpticsPsf(:)), max(psfTargetWFromOpticsOtf(:))]);
    
    % Build up plot
    subplot(6, cols, col);
    imagesc(xSfGridCyclesDegFromOpticsOtf(1, :), ...
        ySfGridCyclesDegFromOpticsOtf(:, 1), ...
        abs(fftshift(otfTargetWFromOptics)), [0 1]);
    axis 'image'
    set(gca, 'XLim', [-60 60], 'YLim', [-60 60]);
    title(sprintf('OTF from optics (%2.0f nm) \n(max = %1.6f)', ...
        targetWavelength, max(otfTargetWFromOptics(:))));
    
    subplot(6, cols, cols + col);
    imagesc(xSfGridCyclesDegFromWvfOtf(1, :), ...
        ySfGridCyclesDegFromWvfOtf(:, 1), ...
        abs(fftshift(otfTargetWFromWvf)), [0 1]);
    axis 'image'
    set(gca, 'XLim', [-60 60], 'YLim', [-60 60]);
    title(sprintf('OTF from wvf (%2.0f nm) \n(max = %1.6f)', ...
        targetWavelength, max(otfTargetWFromWvf(:))));
    
    subplot(6, cols, 2 * cols + col);
    imagesc(xGridMinutesFromWvfPsf(1, :), yGridMinutesFromWvfPsf(:, 1), ...
        psfTargetWFromWvfPsf / maxAll, [0 1]);
    axis 'image'
    set(gca, 'XLim', [-3 3], 'YLim', [-3 3]);
    title(sprintf('PSF (%2.0f nm) direct from wvf\n(sum = %1.6f)', ...
        targetWavelength, sum(psfTargetWFromWvfPsf(:))));
    
    subplot(6, cols, 3 * cols + col);
    imagesc(xGridMinutesFromOpticsPsf(1, :), ...
        yGridMinutesFromOpticsPsf(:, 1), ...
        psfTargetWFromOpticsPsf / maxAll, [0 1]);
    axis 'image'
    set(gca, 'XLim', [-3 3], 'YLim', [-3 3]);
    title(sprintf('PSF (%2.0f nm) direct from optics\n(sum = %1.6f)', ...
        targetWavelength, sum(psfTargetWFromOpticsPsf(:))));
    
   subplot(6, cols, 4 * cols + col);
    imagesc(xGridMinutesFromWvfOtf(1, :), yGridMinutesFromWvfOtf(:, 1), ...
        psfTargetWFromWvfOtf / maxAll, [0 1]);
    axis 'image'
    set(gca, 'XLim', [-3 3], 'YLim', [-3 3]);
    title(sprintf('PSF (%2.0f nm) via wvf OTF\n(sum = %1.6f)', ...
        targetWavelength, sum(psfTargetWFromWvfOtf(:))));
    
    subplot(6, cols, 5 * cols + col);
    imagesc(xGridMinutesFromOpticsOtf(1, :), ...
        yGridMinutesFromOpticsOtf(:, 1), ...
        psfTargetWFromOpticsOtf / maxAll, [0 1]);
    axis 'image'
    set(gca, 'XLim', [-3 3], 'YLim', [-3 3]);
    title(sprintf('PSF (%2.0f nm) via optics OTF\n(sum = %1.6f)', ...
        targetWavelength, sum(psfTargetWFromOpticsOtf(:))));
end
