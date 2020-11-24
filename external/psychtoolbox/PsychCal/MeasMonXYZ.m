function XYZ = MeasMonXYZ(window, settings, whichMeterType)
% XYZ = MeasMonXYZ(window, settings [, whichMeterType=1])
%
% Measure the XYZ of a series of monitor settings.
%
% Usage:
%
% 'window' Onscreen window handle for window to present into.
% 'settings' a 3-rows by nMeas columns matrix with each column defining one
% [r,g,b]' color value to measure.
%
% 'whichMeterType' type of Colorimeter to use:
%
%   0 - Return random spectra. This is useful for testing when you don't have a meter.
% > 0 - Use MeasXYZ(whichMeterType) to measure. See "help MeasXYZ" for available meter types.
%

% History:
%
% 10/26/93  dhb     Wrote it based on ccc code.
% 11/12/93  dhb     Modified to use SetColor.
%  6/23/94  ccc     Modified it from MeasMonSpd.m for
%                   the purpose of measuring XYZ
%  8/9/94   dhb     Added code to go into sync mode
%                   And then commented it out.
% 8/11/94   dhb     Sync mode back in
% 8/15/94   dhb     Sync mode as argument.
% 4/12/97   dhb     New toolbox compatibility, take window and bits args.
% 8/26/97   dhb, pbe Add noMeterAvail option.
% 4/7/99    dhb     Add argument for radius board. Compact default arg code.
% 8/14/00   dhb     Call to CMETER('SetParams') conditional on OS9.
% 8/20/00   dhb     Remove bits arg from call to SetColor.
% 8/21/00   dhb     Remove dependence on RADIUS flag.  This is now handled inside of SetColor.
%           dhb     Change calling conventions to remove unused args.
% 9/14/00   dhb     Sync mode no longer used.  Arg passed for backwards compatibility.
% 2/27/02   dhb, ly Pass whichMeterType rather than noMeterAvail.
% 5/19/20   mk      Fix this disastrous mess somewhat. Still untested.
% 9/02/20   mk      Formatting fixes.

% Check args and make sure window is passed right.
usageStr = 'XYZ = MeasMonXYZ(window, settings [, whichMeterType=1])';
if (nargin < 2 || nargin > 3 || nargout > 1)
    error(usageStr);
end

if ~isscalar(window) || ~isreal(window) || (Screen('WindowKind', window) ~= 1)
    error(usageStr);
end

% Default to PR-650:
if (nargin < 3 || isempty(whichMeterType))
    whichMeterType = 1;
end

[~, nMeas] = size(settings);
XYZ = zeros(3, nMeas);
for i=1:nMeas
    % Set color:
    Screen('FillRect', window, settings(:,i)');
    Screen('Flip', window);

    % Make the measurement:
    switch (whichMeterType)
        case 0
            XYZ(:,i) = sum(settings(:,i))*ones(3,1);
            WaitSecs(0.1);
        otherwise
            XYZ(:,i) = MeasXYZ(whichMeterType);
    end
end
