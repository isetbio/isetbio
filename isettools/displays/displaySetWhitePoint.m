function vci = displaySetWhitePoint(vci, format)
% Set the display white point chromaticity
%
% Syntax:
%   vci = displaySetWhitePoint(vci, format)
%
% Description:
%    The display spectral power distributions of the primaries so that
%    the entered chromaticity is the white point of the display.
%
% Inputs:
%    vci    - Object. The image object
%    format - (Optional) String. Image format. Default 'xyz'.
%
% Outputs:
%    vci    - Object. the modified image object.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: These should be part of displaySet() and displayGet()
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    05/08/18  jnm  Formatting

% Examples:
%{
    % ETTBSkip - vci not initialized
    vci = displaySetWhitePoint(vci, 'xyz')
%}

if ~exist('vci', 'var') || isempty(vci)
    [val, vci] = vcGetSelectedObject('VCIMAGE');
end
if ~exist('format', 'var') || isempty(format), format = 'xyz'; end

switch lower(format)
    case 'xyz'
        wp = chromaticity(imageGet(vci, 'whitepoint'));
        wpchromaticity = ieReadMatrix(wp, '%.3f', ...
            'Enter display whitepoint chromaticity (xy): ');
        if isempty(wpchromaticity), return; end

        Yw = imageGet(vci, 'maxdisplayluminance');
        XYZw = xyy2xyz([wpchromaticity(1), wpchromaticity(2), Yw]);
        displayXYZ = imageGet(vci, 'displayXYZ');

        % XYZw ...
        %   = [1, 1, 1] * diag(sFactor) * displayXYZ ...
        %   = sFactor * displayXYZ
        sFactor = XYZw * inv(displayXYZ);
        spd = imageGet(vci, 'spd');
        vci = vcimageClearData(vci);
        vci.display.spd = spd * diag(sFactor);
        % chromaticity(imageGet(vci, 'whitepoint'))
    otherwise
        error('Unknown format for white point.');
end

return;
