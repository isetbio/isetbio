function oi = oiAdjustIlluminance(oi, newLevel, stat)
% Scale optical image mean illuminance to a new level
%
% Syntax:
%   oi = oiAdjustIlluminance(oi, newLevel, [stat])
%
% Description:
%    Adjust the (mean or max) photon level in the optical structure so that
%    the  mean illuminance level is newLevel, rather than the current
%    level. The fields oi.data.illuminance and oi.data.newLevel are updated
%    as well.
%
%    There are examples contained in the code. To access, type 'edit
%    oiAdjustIlluminance.m' into the Command Window.
%
% Inputs:
%    oi       - Struct. The OI structure.
%    newLevel - Numeric. The desired "stat type" illuminance level as an
%               integer. (Usually mean illuminance level).
%    stat     - (Optional) String. The desired statistic of photon level.
%               Options are 'mean' and 'max' (aka 'peak'). Default 'mean'.
%
% Outputs:
%    oi       - The modified input oi structure.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/01/18  jnm  Formatting
%    06/26/19  JNM  Added second example and removed corresponding TODO.

% Examples:
%{
    OI = oiCreate;
    illu = oiCalculateIlluminance(OI);
    OI = oiSet(OI, 'illuminance', illu);
    OI = oiAdjustIlluminance(OI, 10);
    OI = oiAdjustIlluminance(OI, 100, 'max');
%}
%{
    scene = sceneCreate('uniform');
    scene = sceneSet(scene, 'fov', 15);  % Reasonably large
    scene = sceneAdjustLuminance(scene, 10 ^ -10);

    oi = oiCreate;
    % No lens shading
    optics = oiGet(oi, 'optics');
    optics = opticsSet(optics, 'cos4th', 'off');
    oi = oiSet(oi, 'optics', optics);
    oi = oiCompute(oi, scene);

    [I, meanI, mcI] =  oiCalculateIlluminance(oi);
    OI = oiAdjustIlluminance(oi, 10, 'max');
    I2 = oiCalculateIlluminance(OI);
%}

if notDefined('stat'), stat = 'mean'; end

% Get current OI illuminance
illuminance = oiGet(oi, 'illuminance');

switch lower(stat)
    case {'mean'}
        currentLevel = mean(illuminance(:));
    case {'max', 'peak'}
        currentLevel = max(illuminance(:));
    otherwise
        errordlg('Unknown statistic');
end

s = newLevel / currentLevel;

photons = oiGet(oi, 'photons');
oi = oiSet(oi, 'photons', photons .* s);
oi = oiSet(oi, 'illuminance', illuminance .* s);

end