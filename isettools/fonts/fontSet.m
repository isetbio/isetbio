function font = fontSet(font, param, val, varargin)
% Set the font parameters and derived properties
%
% Syntax:
%   font = fontSet(font, parm, val, [varargin])
%
% Description:
%    Set the font parameters and derived properties
%
% Inputs:
%    font  - The font structure to modify
%    param - The parameter to modify. Possible parameters include:
%            type (Always 'font')
%            name      - String. Default is amalgation of letter, family,
%                        size, and dpi characteristics
%            character - Char. The letter. Default 'g'
%            size      - Integer. Default 14
%            family    - Font family. Default 'Georgia'
%            style     - Font style. Default is 'Normal'.
%            dpi       - Integer. Dots per inch. Default 96
%            bitmap    - Bitmap image. The bitmap of the font is set so the
%                        letter is black on a white background. Use 1 -
%                        bitmap to flip it to white on a black background.
%
% Outputs:
%    font   - The value assigned to the parameter
%
% Optional key/value pairs:
%    **Needs to be filled out**
%
% Notes:
%    * [Note: JNM - Assign someone to fill out optional key/value pairs.]
%

% History:
%    xx/xx/14  HJ/BW  PDCSoft Team, 2014.
%    02/27/18  jnm    Formatting

if notDefined('param'),  error('Parameter must be defined.'); end
if ~exist('val', 'var'), error('Parameter must be defined.'); end

param = ieParamFormat(param);

switch param
    % Book keeping
    case 'type'
        % Always 'font'
    case 'name'
        font.name = val;
    case 'character'
        font.character = val;
    case 'size'
        font.size = val;
    case 'family'
        font.family = val;
    case 'style'
        font.style = val;
    case 'dpi'
        font.dpi = val;
    case 'bitmap'
        font.bitmap = val;
    otherwise
        disp(['Unknown parameter: ', param]);
end

% Always update the font bitmap when you change the size or character.
font.bitmap = fontBitmapGet(font);

end
