function density = coneDensity(eccMM, angle, whichEye, varargin)
% Compute cone packing density as a function of retinal position
%
%   density = coneDensity(ecc, angle, whichEye, varargin)
%
% Inputs:
%   ecc      - eccentricity (retinal position amplitude) in mm
%   angle    - retinal position angle in degree, default is 0
%   whichEye - left or right eye, chosen from 'left' (default) or 'right'
%
% Outputs:
%   density  - cone packing density in cones/mm^2
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
% See also:
%   sensorCreate
%
% HJ, ISETBIO TEAM, 2015

% Check inputs
if notDefined('eccMM'), eccMM = 0; end
if notDefined('angle'), angle = 0; end
if notDefined('whichEye'), whichEye = 'left'; end

% load data
d = load('coneDensity.mat');

% interpolate for retinal position amplitude on axis (nasal, superior,
% temporal and inferior direction)
onAxisD = zeros(5, 1);
angleQ = [0 90 180 270 360];

% compute packing density for superior and inferior
onAxisD(2) = interp1(d.superior.eccMM, d.superior.density, eccMM);
onAxisD(4) = interp1(d.inferior.eccMM, d.inferior.density, eccMM);

% nasal and temporal
switch lower(whichEye)
    case 'left'
        onAxisD(1) = interp1(d.nasal.eccMM, d.nasal.density, eccMM);
        onAxisD(3) = interp1(d.temporal.eccMM, d.temporal.density, eccMM);
    case 'right'
        onAxisD(3) = interp1(d.temporal.eccMM, d.temporal.density, eccMM);
        onAxisD(3) = interp1(d.nasal.eccMM, d.nasal.density, eccMM);
    otherwise
        error('unknown input for whichEye');
end
onAxisD(5) = onAxisD(1);

% Interpolate for angle
density = interp1(angleQ, onAxisD, angle, 'linear');

end