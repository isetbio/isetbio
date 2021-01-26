% Start with a spectral radiance and produce corneal irradiance
%
% Description:
%   Start with spectral radiance and produce corneal irradiance.
%
%   Spectral radiance is specified by a relative spectrum and a luminance
%   in cd/m2.  Corneal radiance comes back as total Watts/m2.  Along the
%   way, spectral corneal irradiance is computed (and then just summed over
%   wavelength).
%
% Require:
%   tbUse('Psychtoolbox-3') or tbUse('isetbio')

% History:
%    12/13/20  dhb  Added some comments.

% Clear
clear; close all;

% Define wavelength spacing.  Use 1 nm spacing
% so that values are per nm.
S = [380 1 401];

% Load D65 spectrum
load spd_D65
theRawSpd = SplineSpd(S_D65,spd_D65,S);

% Load CIE 1931 CMFs
%
% 683 is a magic constant that brings luminance into cd/m2
% when spectral radiance is in Watts/sr-m2-nm
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);

% Specifiy desired lum and scale so that we have spectrum in
% of desired luminance in Watts/sr-m2-nm.
desiredLum = 576;
theSpd = desiredLum*theRawSpd/(T_xyz(2,:)*theRawSpd);

% Compute effective luminous efficiency
theSpd_1 = theRawSpd/sum(theRawSpd);
theLuminance_1 = T_xyz(2,:)*theSpd_1;
fprintf('Luminous efficiency is %d Lumens/Watt\n',theLuminance_1);

% Stimulus size deg
theDegs = 30;
theDegs2 = pi*((theDegs/2)^2);

% Compute corneal irradiance
cornealIrradiance_PowerPerAreaNm = RadianceAndDegrees2ToCornIrradiance(theSpd,theDegs2);

% Total corneal irradiance
cornIrradianceWattsPerArea = sum(cornealIrradiance_PowerPerAreaNm);
fprintf('Corneal irradiance corresponding to %d deg diameter, %d cd/m2 D65: %0.2g Watts/m2\n', ...
        theDegs,desiredLum,cornIrradianceWattsPerArea);
    
% Check by going the other way
checkRadiance = CornIrradianceAndDegrees2ToRadiance(cornealIrradiance_PowerPerAreaNm,theDegs2);
checkLum = T_xyz(2,:)*checkRadiance;
if (abs(desiredLum-checkLum) > 1e-5)
    error('Calculation does not invert as expected');
end
