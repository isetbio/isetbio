function absorptions = computeSingleFrame(obj, oi, varargin)
% Compute function for single frame.
% 
% This function computes mean expected photon absorptions for
% one frame, without including eye movements or noise.
%
% Key/value pairs
%  'fullLMS' - true/false (default false). [WHAT DOES THIS DO?]

% HJ ISETBIO Team 2016

%% Parse inputs
p = inputParser();
p.addRequired('oi', @isstruct);
p.addParameter('fullLMS', false, @islogical);
p.parse(oi, varargin{:});
fullLMS = p.Results.fullLMS;  

%% Get wavelength sampling consistent
%
% Do this by making a copy of current obj and setting wavelength samples to
% be same as oi.
obj = obj.copy();
obj.wave = oiGet(oi,'wave');

%% Get scaled spectral qe
%
% This which includes cone pigment and macular pigment properties. (Lens is
% in oi).  Scale this by the wavelength sample spacing, so that spacing is
% taken into account when we compute isomerizations (aka absorptions in the
% isetbio world.)
sQE = obj.qe * oiGet(oi, 'bin width');

%% Compute cone isomerization density oi sampled locations for each class of cones
%
% These need to be scaled by cone integration area and time to get actual
% absorptions.
%
% Also, there are not necessarily cones at all of these locations, we'll
% sample below.
[photons, r, c] = RGB2XWFormat(oiGet(oi, 'photons'));
absDensityLMS = XW2RGBFormat(photons * sQE, r, c);

% Regrid the isomerizatoin density from oi sample locations to cone locations
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

%% Sometimes we don't have the cone type so we have a bad
% number. Set the missing values to 0
absDensity(isnan(absDensity)) = 0;

% compute expected cone absorptions
absorptions=absDensity*obj.pigment.pdArea*obj.integrationTime;

end