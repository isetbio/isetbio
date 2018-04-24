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
%    newLevel - Integer. The desired "stat type" illuminance level.
%               (Usually mean illuminance level).
%    stat     - (Optional) String. The desired statistic of photon level.
%               Options are 'mean' and 'max' (aka 'peak'). Default 'mean'.
%
% Outputs:
%    oi       - The modified input oi structure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - TODO: Assign someone to fix the example. I've fleshed
%      it out enough to stop the error message, but it does not create a
%      real illuminance, or modify it.]
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/01/18  jnm  Formatting

% Examples:
%{
    OI = oiCreate;
    illu = oiCalculateIlluminance(OI);
   OI = oiSet(OI, 'illuminance', illu);
    OI = oiAdjustIlluminance(OI, 10);
   OI = oiAdjustIlluminance(OI, 100, 'max');
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