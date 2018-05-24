function vci = displaySetMaxLuminance(vci)
% Set simulated display SPD to match a peak luminance
%
% Syntax:
%   vci = displaySetMaxLuminance(vci)
%
% Description:
%    The peak luminance is set separately. Perhaps we should find a way
%    to ensure they can never be inconsistent.
%
% Inputs:
%    vci - Object. The image object.
%
% Outputs:
%    vci - Object. The modified image object.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: These should be part of displaySet() and displayGet()
%

% Examples:
%{
	% ETTBSkip - skipping for imageGet reference (deprecated?)
    vci = displaySetMaxLuminance;
    vci = displaySetMaxLuminance(vci)
%}

if notDefined('vci'), [~, vci] = vcGetSelectedObject('VCIMAGE'); end

Yw = imageGet(vci, 'maxdisplayluminance');

maxLum = ieReadNumber('Enter desired max display luminance (Y): ', ...
    Yw, '%.2f');
if isempty(maxLum), return; end

sFactor = maxLum / Yw;
spd = imageGet(vci, 'spd');
vci = vcimageClearData(vci);
vci = imageSet(vci, 'spd', spd * sFactor);

return;
