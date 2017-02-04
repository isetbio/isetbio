function [xGridMinutes,yGridMinutes,psf] = OtfToPsf(xSfGridCyclesDeg,ySfGridCyclesDeg,otf,varargin)
%OTFTOPSF  Convert a 2D optical transfer fucntion to a 2D point spread function.
%    [xSfGridCyclesDeg,ySfGridCyclesDeg,otf] = PsfToOtf(xGridMinutes,yGridMinutes,psf,varargin)
%
%    Converts a optical transfer function specified over two-dimensional
%    spatial frequency in cycles per degree to a point spread function
%    specified over positions in minutes.  For human vision, these are each
%    natural units.
%
%    The input spatial frequencies should be specified in matlab's grid
%    matrix format and sf for x and y should be specified over the same
%    spatialfrequency extent and with the same number of evenly spaced
%    samples.  The number of samples should be even. Position (0,0) should
%    be at location n/2+1 in each dimension.
%
%    Positions are returned using the same conventions.
%
%    No normalization is performed.  The psf should be real, and we
%    complain (throw an error) if it is not, to reasonable numerial
%    precision. 
%
%    We wrote this rather than simply relying on Matlab otf2psf because
%    we want a routine that deals with the conversion of spatial frequency to
%    spatial support.
%
%    PSYCHOPTICSTEST shows that this works very well when we go back and
%    forth for diffraction limited OTF/PSF.  But not exactly exactly
%    perfectly.  A signal processing maven might be able to track down
%    whether this is just a numerical thing or whether some is some small
%    error, for example in how position is converted to sf or back again in
%    the PSFTOOTF.
%
%    See also PSFTOOTF, PSYCHOPTICSTEST.

%% Reality checks on passed input
[m,n] = size(xSfGridCyclesDeg);
centerPosition = n/2 + 1;
if (m ~= n)
    error('psf must be passed on a square array');
end
if (mod(m,2) ~= 0)
    error('Number of rows of spatial support must be even');
end
if (mod(n,2) ~= 0)
    error('Number of cols of spatial support must be even');
end
[m1,n1] = size(ySfGridCyclesDeg);
if (m1 ~= m || n1 ~= n)
    error('x and y positions are not consistent');
end
[m2,n2] = size(otf);
if (m2 ~= m || n2 ~= n)
    error('x and y positions are not consistent');
end
if (~all(xSfGridCyclesDeg(:,centerPosition) == 0))
    error('Zero spatial frequency is not in right place in the passed xGrid');
end
if (~all(ySfGridCyclesDeg(centerPosition,:) == 0))
    error('Zero spatial frequency is not in right place in the passed yGrid');
end
if (xSfGridCyclesDeg(1,centerPosition) ~= ySfGridCyclesDeg(centerPosition,1))
    error('Spatial frequency extent of x and y grids does not match');
end
diffX = diff(xSfGridCyclesDeg(:,centerPosition));
if (any(diffX ~= diffX(1)))
    error('X positions not evenly spaced');
end
diffY = diff(ySfGridCyclesDeg(centerPosition,:));
if (any(diffY ~= diffY(1)))
    error('Y positions not evenly spaced');
end
if (diffX(1) ~= diffY(1))
    error('Spatial frequency sampling in x and y not matched');e
end

%% Generate position grids
%
% Samples are evenly spaced and the same for both x and y (checked above).
spatialFrequencyExtentCyclesDeg = xSfGridCyclesDeg(centerPosition,end)-xSfGridCyclesDeg(centerPosition,1);
spatialFrequencyExtentCyclesMinute = spatialFrequencyExtentCyclesDeg/60;
positionMinutes = (-n/2:n2/2-1)/spatialFrequencyExtentCyclesMinute;
[xGridMinutes,yGridMinutes] = meshgrid(positionMinutes,positionMinutes);
    
%% Compute otf
psf = fftshift(ifft2(ifftshift(otf)));

%% See if there is stray imaginary stuff, get rid of it if so.
%
% Throw an error if the returned psf isn't in essence real valued.
if (any(abs(imag(psf(:))) > 1e-10))
    error('Computed psf is not sufficiently real');
end
psf = abs(psf);

end