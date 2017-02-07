function [xSfGridCyclesDeg,ySfGridCyclesDeg,otf] = PsfToOtf(xGridMinutes,yGridMinutes,psf,varargin)
%PSFTOOTF  Convert a 2D point spread function to a 2D optical transfer fucntion.
%    [xSfGridCyclesDeg,ySfGridCyclesDeg,otf] = PsfToOtf(xGridMinutes,yGridMinutes,psf,varargin)
%
%    Converts a point spread function specified over two-dimensional
%    positions in minutes to a optical transfer function specified over
%    spatial frequency in cycles per degree.  For human vision, these are
%    each natural units.
%
%    The input positions should be specified in matlab's grid matrix format
%    and x and y should be specified over the same spatial extent and with
%    the same number of evenly spaced samples. Position (0,0) should be at
%    location floor(n/2)+1 in each dimension.
%
%    Spatial frequencies are returned using the same conventions.
%
%    No normalization is performed.  If the phase of the OTF are very small
%    (less than 1e-10) the routine assumes that the input psf was spatially
%    symmetric around the origin and takes the absolute value of the
%    computed otf so that the returned otf is real.
%
%    We wrote this rather than simply relying on Matlab's psf2otf because
%    we don't understand quite how that shifts position of the passed psf
%    and because we want a routine that deals with the conversion of
%    spatial support to spatial frequency support.
%
%    PSYCHOPTICSTEST shows that this works very well when we go back and
%    forth for diffraction limited OTF/PSF.  But not exactly exactly
%    perfectly.  A signal processing maven might be able to track down
%    whether this is just a numerical thing or whether some is some small
%    error, for example in how position is converted to sf or back again in
%    the OTFTOPSF.
%
%    See also OTFTOPSF, PSYCHOPTICSTEST.

%% Reality checks on passed input
[m,n] = size(xGridMinutes);
centerPosition = floor(n/2) + 1;
if (m ~= n)
    error('psf must be passed on a square array');
end
[m1,n1] = size(yGridMinutes);
if (m1 ~= m || n1 ~= n)
    error('x and y positions are not consistent');
end
[m2,n2] = size(psf);
if (m2 ~= m || n2 ~= n)
    error('x and y positions are not consistent');
end
if (~all(xGridMinutes(:,centerPosition) == 0))
    error('Zero position is not in right place in the passed xGrid');
end
if (~all(yGridMinutes(centerPosition,:) == 0))
    error('Zero position is not in right place in the passed yGrid');
end
if (xGridMinutes(1,centerPosition) ~= yGridMinutes(centerPosition,1))
    error('Spatial extent of x and y grids does not match');
end
diffX = diff(xGridMinutes(:,centerPosition));
if (any(diffX ~= diffX(1)))
    error('X positions not evenly spaced');
end
diffY = diff(yGridMinutes(centerPosition,:));
if (any(diffY ~= diffY(1)))
    error('Y positions not evenly spaced');
end
if (diffX(1) ~= diffY(1))
    error('Spatial sampling in x and y not matched');e
end

%% Generate spatial frequency grids
[xSfGridCyclesDeg,ySfGridCyclesDeg] = PositionGridMinutesToSfGridCyclesDeg(xGridMinutes,yGridMinutes);
%% Compute otf
otf = fftshift(fft2(ifftshift(psf)));

%% See if there is stray imaginary stuff, get rid of it if so.
%
% We don't require that the input psf be symmetric, so there could be
% actual imaginary values.  Thus we do our best to make a good guess.
if (all(abs(imag(otf(:))) < 1e-10))
    otf = abs(otf);
end

end