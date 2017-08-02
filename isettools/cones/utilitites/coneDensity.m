function density = coneDensity(ecc, angleDeg, whichEye, varargin)
%coneDensity  Compute cone packing density as a function of retinal position
%
% Usage:
%     density = coneDensity(ecc, angle, whichEye, varargin)
%
% Description:
%     Compute cone packing density as a function of retinal position.
%
% Inputs:
%     ecc      - retinal eccentricity in meters, default is 0 (.30 mm/deg)
%     angle    - angle in degree, default is 0
%     whichEye - left or right eye, chosen from 'left' (default) or 'right'
%
% [DHB NOTE: WHAT IS ANGULAR COORDINATE SYSTEM? WHICH WAY IS ZERO DEGS, AND 
% WHICH WAY DOES ANGLE GO.]
%
% Outputs:
%     density  - cone packing density in cones/mm^2
%
% References:
%   1) Curcio, C. A., Sloan, K. R., Kalina, R. E. and Hendrickson, A. E.
%      (1990), Human photoreceptor topography. J. Comp. Neurol., 292: 
%      497?523. doi: 10.1002/cne.902920402
%   2) Song, H., Chui, T. Y. P., Zhong, Z., Elsner, A. E., & Burns, S. A.
%      (2011). Variation of Cone Photoreceptor Packing Density with Retinal
%      Eccentricity and Age. Investigative Ophthalmology & Visual Science,
%      52(10), 7376?7384. http://doi.org/10.1167/iovs.11-7199
%
% Example:
%   density = coneDensity(8*1e-3, 10, 'left');
%
% See also:
%   coneMosaic
%
% HJ, ISETBIO TEAM, 2015

% Check inputs
if notDefined('ecc'), ecc = 0; end
if notDefined('angle'), angleDeg = 0; end
if notDefined('whichEye'), whichEye = 'left'; end

% load cone density data.  Units are cones / mm^2
d = load('coneDensity.mat');

% interpolate for retinal position amplitude on axis (nasal, superior,
% temporal and inferior direction)
onAxisD = zeros(5, numel(ecc));
angleQ = [0 90 180 270 360];

% Convert to mm for functions below
eccMM = ecc*1e3;

% compute packing density for superior and inferior
onAxisD(2,:) = interp1(d.superior.eccMM, d.superior.density, eccMM);
onAxisD(4,:) = interp1(d.inferior.eccMM, d.inferior.density, eccMM);

% nasal and temporal
switch lower(whichEye)
    case 'left'
        onAxisD(1,:) = interp1(d.nasal.eccMM, d.nasal.density, eccMM);
        onAxisD(3,:) = interp1(d.temporal.eccMM, d.temporal.density, eccMM);
    case 'right'
        onAxisD(1,:) = interp1(d.temporal.eccMM, d.temporal.density, eccMM);
        onAxisD(3,:) = interp1(d.nasal.eccMM, d.nasal.density, eccMM);
    otherwise
        error('unknown input for whichEye');
end
onAxisD(5,:) = onAxisD(1,:);

% Interpolate for angle
density = interp1(angleQ, onAxisD, angleDeg, 'linear');

end
