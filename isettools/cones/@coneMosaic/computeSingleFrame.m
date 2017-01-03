function absorptions = computeSingleFrame(obj, oi, varargin)
% computeSingleFrame  Single frame compute function for coneMosaic object.
%     absorptions = computeSingleFrame(obj, oi, varargin)
%
% This function computes mean expected photon isomerizations (also called
% absorptions in isetbio) for one frame, without including eye movements or
% noise.
%
% Key/value pairs
%  'fullLMS' - true/false (default false). Return values for a full mosaic,
%    that is for mosaic with L, M, and S cones at each cone position.  This
%    is row by col by 3 matrix, where row and column are the mosaic
%    dimensions.  This is not biologically realistic but useful for some
%    computations.
  
% HJ ISETBIO Team 2016

%% Parse inputs
p = inputParser();
p.addRequired('oi',@isstruct);
p.addParameter('fullLMS', false,@islogical);
p.parse(oi, varargin{:});
fullLMS = p.Results.fullLMS;  

%% Get wavelength sampling consistent
%
% Do this by making a copy of current obj and setting wavelength samples to
% be same as oi.
obj = obj.copy();
obj.wave = oiGet(oi,'wave');

%% Get wavelength-spacing scaled spectral qe
%
% This which includes cone pigment and macular pigment properties. (Lens is
% in oi).  Scale this by the wavelength sample spacing, so that spacing is
% taken into account when we compute isomerizations (aka absorptions in the
% isetbio world.)
sQE = obj.qe * oiGet(oi,'bin width');

%% Compute cone isomerization density oi sampled locations for each class of cones
%
% These need to be scaled by cone integration area and time to get actual
% isomerizations.
%
% Note that there are not necessarily cones at all of these locations, we'll
% sample below.
[photons, r, c] = RGB2XWFormat(oiGet(oi, 'photons'));
absorbDensityLMS = XW2RGBFormat(photons * sQE, r, c);

%% Blur by cone aperture
%
% Physically, this operation occurs as the light is gathered by the cone,
% before integrating over wavelength.  We treat the cone aperature as
% independent of wavelength, however, so the convolution commutes with the
% summation over wavelength used above to get the isomerization density for
% each class of cone.  It's faster to convolve here, since there are fewer
% bands to deal with.
if (obj.apertureBlur)
    % Make the blur kernal.
    %
    % First convert area of cone aperture to a radius
    apertureRadius = sqrt(obj.pigment.pdArea/pi);
    
    % Get optical image resolution and figure out how big conv kernal needs
    % to be.  Make sure it is an odd number greater than 0.
    oiRes = oiGet(oi, 'height spatial resolution');
    if (oiRes ~= oiGet(oi, 'width spatial resolution'))
        error('Cannot do blurring by cone apearture if oi row and column resolutions differ');
    end
    apertureSamples = ceil(2*apertureRadius/oiRes);
    if (mod(apertureSamples,2) == 0)
        apertureSamples = apertureSamples+1;
    end
    
    % Make the kernal and make the circular pillbox with unit volume.
    apertureKernal = zeros(apertureSamples,apertureSamples);
    [apertureRows,apertureCols] = meshgrid(sample2space(1:apertureSamples,1:apertureSamples,oiRes,oiRes));
    apertureRadii = sqrt(apertureRows.^2+apertureCols.^2);
    apertureKernal(apertureRadii <= apertureRadius) = 1;
    apertureKernal = apertureKernal./sum(apertureKernal(:));
 
    % Do the convolution for each cone type.
    %
    % In the loop over cone types, 1 means blank/black so we just iterate 2:4.
    for ii = 2:4
        absorbDensityLMS(:,:,ii-1) = conv2(absorbDensityLMS(:,:,ii-1),apertureKernal,'same');
    end 
end

% Regrid the isomerization density from oi sample locations to cone
% locations.  The optical image and mosaic are assumed to be centered on
% one another, with all positions treated relative to the center.
%
% oiR and oiC give row and column positions of each spatial sample in the
% oi in microns relative to the center of the oi.
%
% The original args to this were 0:r-1 and 0:c-1, but the way that
% sample2space is written, it doesn't make any difference and that seemed 
% more confusing.
[oiR, oiC] = sample2space(1:r, 1:c, ...
    oiGet(oi, 'height spatial resolution'), ...
    oiGet(oi, 'width spatial resolution'));

% coneR and coneC give the row and column positions of each cone in microns
% relative to the center of the mosaic.
%
% The original args to this were 0.5:r-0.5 and 0.5:r-0.5, but the way that
% sample2space is written, it doesn't make any difference and that seemed 
% really really confusing because it gave the impression that the cone
% mosaic was treated differently from the optical image above.
[coneR, coneC] = sample2space(1:obj.rows,1:obj.cols, ...
    obj.patternSampleSize(2), obj.patternSampleSize(1));

% Allocate space for computed density.  If we are returning the LMS values
% at each cone position, independent of type, it's row by col by 3.
if fullLMS
    absorbDensity = zeros(obj.rows,obj.cols,3);
else
    absorbDensity = zeros(obj.rows,obj.cols);
end

% Loop through L, M and S cones and get isomerizations.
%
% In the loop over cone types, 1 means blank/black so we just iterate 2:4.
warning('off','MATLAB:interp1:NaNinY');
for ii = 2:4
    % Get the isomerization density for the current cone class at the row
    % positions of all the cones.  Note that interp1 operates on the first
    % dimension of the matrix passed as the second argument, independent of
    % whether the first and third entries are row or column vectors.
    absorbDensityOneConeClass = interp1(oiR, absorbDensityLMS(:,:,ii-1), ...
        coneR, 'linear', 0)';
    
    % Now interpolate in place across rows to get the values at the column
    % positions, so that we have the density at each cone position.
    absorbDensityOneConeClass = interp1(oiC, absorbDensityOneConeClass, coneC, 'linear', 0)';
    
    % Save density for each cone class.  Either the whole array of we're
    % returning the full LMS isomerizations, or one value for each cone
    % location if we are in the usual case of simulating the interleaved
    % mosaic.
    if fullLMS
        absorbDensity(:,:,ii-1) = absorbDensityOneConeClass;
    else
        absorbDensity = absorbDensity+(obj.pattern==ii).*absorbDensityOneConeClass;
    end
end
warning('on','MATLAB:interp1:NaNinY');

% Sometimes we don't have the cone type at location so we have a bad number. Set the missing values to 0
absorbDensity(isnan(absorbDensity)) = 0;

% Multiply by integration area and time to convert isomerization density to
% isomerizations.
absorptions=absorbDensity*obj.pigment.pdArea*obj.integrationTime;

end