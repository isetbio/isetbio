function OImax = oiExtractBright(oi)
% Extract brightest pixel from an optical image and return it as an OI
%
% Syntax:
%   OImax = oiExtractBright([oi])
%
% Description:
%    Make a 1 pixel optical image with a spectral illumination of the
%    highest illuminance pixel. The optics settings are adjusted so that
%    OTF and off-axis computations are skipped. This is used in setting
%    exposure duration (autoExposure) and other OI evaluations.
%
%    There are examples listed in the code below. To access, type 'edit
%    oiExtractBright.m' into the Command Window.
%
% Inputs:
%    oi    - (Optional) Struct. An optical image structure. Default is to
%            pull in an existing OI.
%
% Outputs:
%    OImax - Struct. An optical image structure containing the brightest
%            pixel from the provided oi.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - TODO: Assign someone to fix example. Irradiance has to
%      be calculated for illuminance to calculate, and then be assigned
%      back. Currently I've only added the required scene and oi that are
%      needed to calculate a base irradiance and luminance. Current error
%      is "Index exceeds matrix dimensions. Error in oiCalculateIlluminance
%      (line 100) irradianceP = oiGet(oi, 'photons', wave(ii));"]
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/06/18  jnm  Formatting

% Example:
%{
    % ETTBSkip - Example erroring with 'exceeds matrix dimensions'
    scene = sceneCreate;
    oi = oiCreate;
    oi = oiSet(oi, 'photons', oiCalculateIrradiance(scene, oi))
    oi = oiSet(oi, 'illuminance', oiCalculateIlluminance(oi))
    OImax = oiExtractBright(oi);
%}

if notDefined('oi'), oi = vcGetObject('OI'); end

% Find the brightest part of the scene
sz = oiGet(oi, 'size');
illuminance = oiGet(oi, 'illuminance');

% If illuminance has not been computed, compute it here
[~, ind] = max(illuminance(:));
[rect(2), rect(1)] = ind2sub(sz, ind);
rect(3) = 1;
rect(4) = 1;

% Now, we crop the data to form a small opticalimage containing only the
% highest illuminance.
OImax = oiCrop(oi, rect);

% We adjust the optics
optics = oiGet(OImax, 'optics');
optics = opticsSet(optics, 'otfmethod', 'skip');
optics = opticsSet(optics, 'offaxismethod', 'skip');

OImax = oiSet(OImax, 'optics', optics);

end