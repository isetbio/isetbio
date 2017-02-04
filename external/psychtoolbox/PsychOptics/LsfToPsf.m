function psf = LsfToPsf(lsf,varargin)
%LSFToPSF  Convert a line spread function to a point spread function.
%   psf = LSFTOPSF(lsf,varargin)
%
%   This works by taking the lsf into the one-dimensional frequency domain
%   (i.e. get the 1D MTF) and then creating a circularly symmetric version
%   of this.  The 2D frequency representation is then converted back to the
%   spatial domain to produce the psf.  This produces a spatially-symmetric
%   PSF consistent with the measured line spread function.
%
%   This method is described by Marchand, 1964, JOSA, 54, 7, pp. 915-919
%   and is one of several methods provided.  In 1964, taking the Fourier
%   transform was computationally intense. 
%
%   There is a second paper by Marchand (1965, JOSA, 55, 4, 352-354) which
%   treats the more general case where you have line spread functions for
%   many orientations and want to recover an psf that is not necessarily
%   spatially-symmetric.
%
%   The lsf must be spatially symmetric.  This makes sense given that we
%   are going to recover a spatially symmetric psf.
%
%   The lsf must be defined on a support with an even number of samples
%   that have the general feature that the zero position is at location
%   n/2 + 1, where n is the number of samples. This means that there is one
%   more negative position than positive position in the spatial support.
%   Doing it this way is important because of the way the FFT works.
%   Probably you could write code to do it with odd support, but that would
%   require more fussing than I am willing to do.
%
%   You want to make sure that the spatial suport is large enough to
%   capture the full lsf.
%
%   The returned psf is normalized to have unit volume.
%
%   See also PSFTOLSF, PSYCHOPTICSTEST

% Check for even support and symmetry
n1D = length(lsf);
if (mod(n1D,2) ~= 0)
    error('LSF must be defined on spatial support with even dimension.');
end
centerPosition = n1D/2 + 1;
negValues = lsf(centerPosition-1:-1:2);
posValues = lsf(centerPosition+1:end);
if (any(posValues ~= negValues))
    error('LSF is not spatially symmetric');
end
    
% Convert LSF to 1-D MTF
%
% Taking the absolute value.  The fftshift calls wrapping the fft lets us
% work in a coordinate system where 0 in both space and frequency is at
% position n1D/2+1.  If n1D were odd, we'd have to figure out more
% precisely the difference between fftshift and ifftshift, this may not be
% used completely correctly here.  My current view is that ifftshift takes
% one from 0 centered to 0 at first entry, and fftshift goes the other way.
lsfNormalized = lsf/sum(lsf(:));
lsfMTFCenteredRaw = fftshift(fft(ifftshift(lsf)));

% Result should be real.  Check that it's OK to tolerance and then get rid
% of any pesky imaginary components.
lsfOTFCenteredImag = imag(lsfMTFCenteredRaw);
if (any(abs(lsfOTFCenteredImag > 1e-10)))
    error('Imaginary part of 1D MTF is bigger than makes sense from numerical error alone');
end
lsfOTFCentered = abs(lsfMTFCenteredRaw);

% Create symmetric 2D MTF from 1D MTF
%
% Values outside the frequency support of the lsf-based OTF are set to
% zero.  This should be fine as long as the spatial support is large enough
% so that the MTF has fallen to zero at highest frequencies specified.
radiusMatrix = MakeRadiusMat(n1D,n1D,centerPosition);
psfOTFCentered = interp1(0:n1D/2-1,lsfOTFCentered(centerPosition:end),radiusMatrix,'linear',0);

% Do inverse fft, now in 2D.  Again wrap with ifftshift and fftshift so that we
% are working in 0 centered coordinates.
psfRaw = fftshift(ifft2(ifftshift(psfOTFCentered)));

% Result should be real.  Check that it's OK to tolerance and then get rid
% of any pesky imaginary components.
psfImag = imag(psfRaw);
if (any(abs(psfImag) > 1e-10))
    error('Imaginary component of derived psf too big for comfort.');
end
psf = abs(psfRaw);

% Normalize to unit volume
psf = psf/sum(psf(:));

end

