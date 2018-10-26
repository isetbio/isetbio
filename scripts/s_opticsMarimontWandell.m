% s_OpticsMarimontWandell
%
% Description:
%    This script shows the Marimont-Wandell OTF and corresponding PSF.
%    It also explores the difference between the original way that
%    isetbio computed the psf and the current method.
%
%    The current method follows Matlab FFT conventions and seems like the
%    right thing to do, and the two methods give similar, but non-the-less
%    noticably different, results for the PSF.
%
%    The plot of the OTF looks quite similar to Figure 3 of M-W.
%
%    A plot of the LSF is not so close to Figure 4 of the M-W paper, with a
%    "peakier" LSF obtained by the current code. The LSF is computed from
%    the OTF using code separate from that which computes the PSF.  The LSF
%    computation is in oiPlot, and looks like it uses the right Matlab
%    conventions.
%
%    The negative values in the original "incorrect" method are smaller
%    than those in the current "correct" method. I am not sure if that has
%    any significance.
%
% 03/31/18  dhb  Wrote it.
% 10/26/18  dhb  Another try.

%% Clear
clear; close all;

%% Get Marimont-Wandell optics
%
% These are computed as the OTF, following
% the M-W paper.  Routine humanCore is key,
% see also humanAchromaticOTF.  The OTF has
% 1 (DC) in the upper left at each wavelength,
% as we expect.
oi = oiCreate('human');
optics = oiGet(oi,'optics');
OTFMarimontWandell = optics.OTF;

%% Plot the OTF
%
% Compare with Figure 3 of M-W. Note that
% here frequency is in cycles per mm and that
% 1 cyc/deg is about 3 cycles mm.  Set the
% axis cutoff to 100 cycles/mm to correspond
% to the 30 cycles/deg in Figure 3.
oiPlot(oi,'otf wavelength');
xlim([0 100]);
view([31.5 25.5]);

%% Plot the LSF
%
% Compare with Figure 4 of M-W. 
% Axis of +/- 150 um is about +/-
% 0.5 deg.
oiPlot(oi,'ls wavelength');
xlim([-150 150]);
zlim([-0.1 0.6]);
view([19.5 19.5]);

% Get number of samples
nSamples = size(OTFMarimontWandell.OTF,1);
centerSample = floor(nSamples/2)+1;

%% Convert to PSF the ISETBIO way
% 
% This now uses OtfToPsf to do the conversion, with 
% fftshift inserted to put DC in the center before the
% call. There are some negative
% values in the PSF, so the call to OtfToPsf is done
% with a sufficiently large tolerance not to throw an error.
%
% Make sure it is normalized to unit volume at each wavelength.
PSFISETBio = opticsGet(optics,'psf data');
maxFracNegative = 0;
for ww = 1:size(PSFISETBio,3)
    volume = sum(sum(PSFISETBio(:,:,ww)));
    if (abs(volume-1) > 1e-10)
        error('PSF does not have volume of 1');
    end
    PSFISETBio(:,:,ww) = PSFISETBio(:,:,ww)/volume;
    minPSF1Val = min(min(PSFISETBio(:,:,ww)));
    if (minPSF1Val < -1e-15)
        fprintf('PSF1 has negative value for wl index %d, value is %g\n',ww,minVal);
    end
    
    % This is what we are actually doing at each wavelength in the above,
    % see OtfToPsf.  Note that the raw psf has unit volume, which it
    % should.
    temp = OTFMarimontWandell.OTF(:,:,ww);
    temp1 = ifftshift(fftshift(temp));
    if (any(temp ~= temp1))
        error('ifftshift does not invert fftshift');
    end
    psf = fftshift(ifft2(ifftshift(fftshift(temp))));
    if (abs(sum(psf(:))-1) > 1e-10)
        fprintf('Raw ISETBio PSF does not have unit volume\n');
    end
    if (any(abs(imag(psf(:))) > 1e-10))
        error('Computed psf is not sufficiently real');
    end
    if (any(imag(psf(:))) ~= 0)
        psf = psf - imag(psf)*1i;
    end
    if (max(psf(:)) <= 0)
        error('Computed PSF has no positive values.  This is not good.');
    end
    if (min(psf(:)) < 0 && abs(min(psf(:)))/max(psf(:)) > maxFracNegative)
        maxFracNegative = abs(min(psf(:)))/max(psf(:));
    end
    psf(psf < 0) = 0;
    psf = abs(psf);
    psf = psf/sum(psf(:));
    
    % Check that we get the same answer
    psfISETBio = PSFISETBio(:,:,ww);
    if (max(abs(psfISETBio(:) - psf(:))) > 1e-6*max(psfISETBio(:)))
        error('Cannot do same calculation two ways');
    end
    PSFCheck(:,:,ww) = psf;

end
fprintf('Maximum fractional negative PSF in ISETBio method: %f\n',maxFracNegative);

%% Do it manually, the way it was coded.
%
% The original isetbio implementation used fft2 rather than ifft2.
% And the OTF was stored from DC at center using fftshift rather than
% ifftshift.  The two fftshifts in the arg to fft below reproduce the
% earlier fftshift behavior.
maxFracNegative = 0;
for ww = 1:size(OTFMarimontWandell.OTF,3)
    % Get PSF the original isetbio way. Note that the raw psf does not have
    % unit volumen with this method.
    PSFOrig(:,:,ww) = fftshift(fft2(fftshift(fftshift(OTFMarimontWandell.OTF(:,:,ww)))));
    if (abs(sum(sum(PSFOrig(:,:,ww)))-1) < 1e-3)
        fprintf('Orig PSF does have unit volume\n');
    end
    if (min(min(PSFOrig(:,:,ww))) < 0 & abs(min(min(PSFOrig(:,:,ww)))/max(max(PSFOrig(:,:,ww)))) > maxFracNegative)
        maxFracNegative = abs(min(min(PSFOrig(:,:,ww)))/max(max(PSFOrig(:,:,ww))));
    end
    
    % Check for imag values and zero out imag component.
    % This was not done in the original isetbio implementation
    if (any(abs(imag(PSFOrig(:,:,ww))) > 1e-10))
        error('Computed psf is not sufficiently real');
    end
    if (any(imag(PSFOrig(:,:,ww))) ~= 0)
        PSFOrig(:,:,ww) = PSFOrig(:,:,ww) - imag(PSFOrig(:,:,ww))*1i;
    end
    
    % Remove negative values and normalize. Could worry here about
    % imaginary values as well.
    temp = PSFOrig(:,:,ww);
    temp(temp < 0) = 0;
    PSFOrig(:,:,ww) = abs(PSFOrig(:,:,ww));
    PSFOrig(:,:,ww) = PSFOrig(:,:,ww)/sum(sum(PSFOrig(:,:,ww)));
end
fprintf('Maximum fractional negative PSF in original method: %f\n',maxFracNegative);

%% Check for agreement at each between new and old ways of getting PSF
%
% These bear reasonable but not perfect resemblence to each other. The old
% way does not generate negative PSF values, whereas the new way does.
% This is a little mysterious to me.
maxFracMismatch = 0;
maxMismatchWlIndex = 0;
for ww = 1:size(OTFMarimontWandell.OTF,3)
    PSF1ISETEBioThisWl = PSFISETBio(:,:,ww);
    PSFOrigThisWl = PSFOrig(:,:,ww);
    fracMismatch(ww) = max(abs(PSFOrigThisWl(:)-PSF1ISETEBioThisWl(:)))/max(PSF1ISETEBioThisWl(:));
    if (fracMismatch(ww) > maxFracMismatch)
        maxFracMismatch = fracMismatch(ww);
        maxMismatchWlIndex = ww;
    end
    
    % Plot a slice through PSF
    figure; clf; hold on
    plot(PSF1ISETEBioThisWl(centerSample,:),'r','LineWidth',3);
    plot(PSFOrigThisWl(centerSample,:),'b','LineWidth',2);
    legend({'ISETBio PSF slice','Old method PSF slice'},'Location','NorthEast');
    xlabel('Position');
    ylabel('PSF');
    title(sprintf('Wavelength %d',OTFMarimontWandell.wave(ww)));
end
drawnow;
fprintf('Maximum fractional mismatch between current and original method: %f at wl index %d\n',maxFracMismatch,maxMismatchWlIndex);
[~,index] = sort(fracMismatch,'ascend');
for ww = 1:size(OTFMarimontWandell.OTF,3)
    figure(index(ww));
end












