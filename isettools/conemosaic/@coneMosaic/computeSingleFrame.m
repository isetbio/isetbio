function absorptions = computeSingleFrame(obj, oi, varargin)
% Compute function for single frame.
% 
% This function computes mean expected photon absorptions for
% one frame, without including eye movements or noise.
% 
% HJ ISETBIO Team 2016

% parse inputs
p = inputParser();
p.addRequired('oi', @isstruct);
p.addParameter('fullLMS', false, @islogical); % full LMS images

p.parse(oi, varargin{:});

fullLMS = p.Results.fullLMS;  % Logical

% make a copy of current obj and set cone wavelength samples to
% be same as oi
obj = obj.copy();
obj.wave = oiGet(oi, 'wave');

% get scaled spectral qe, which includes cone pigment and
% macular pigment properties. (Lens is in oi).
sQE = obj.qe * oiGet(oi, 'bin width');

% compute cone absorption density at oi sampled locations
[photons, r, c] = RGB2XWFormat(oiGet(oi, 'photons'));
absDensityLMS = XW2RGBFormat(photons * sQE, r, c);

% regrid the density from oi sample locations to cone locations
[oiR, oiC] = sample2space(0:r-1, 0:c-1, ...
    oiGet(oi, 'height spatial resolution'), ...
    oiGet(oi, 'width spatial resolution'));
[coneR, coneC] = sample2space(0.5:obj.rows - 0.5, ...
    0.5:obj.cols - 0.5, obj.patternSampleSize(2), obj.patternSampleSize(1));

if fullLMS
    absDensity = zeros(obj.rows, obj.cols, 3);
else
    absDensity = 0;
end
warning('off','MATLAB:interp1:NaNinY');
for ii = 2 : 4  % loop through L, M and S, 1 = Blank/Black
    curDensity = interp1(oiR, absDensityLMS(:,:,ii-1), ...
        coneR, 'linear', 0)';
    curDensity = interp1(oiC, curDensity, coneC, 'linear', 0)';
    if fullLMS
        absDensity(:, :, ii-1) = curDensity;
    else
        % Pick out the relevant cones by position
        absDensity = absDensity+(obj.pattern==ii).*curDensity;
    end
end
warning('on','MATLAB:interp1:NaNinY');

% Sometimes we don't have the cone type so we have a bad
% number. Set the missing values to 0
absDensity(isnan(absDensity)) = 0;

% compute expected cone absorptions
absorptions=absDensity*obj.pigment.pdArea*obj.integrationTime;
end