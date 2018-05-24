function val = fontGet(font, param, varargin)
% Get font parameters and derived properties 
%
% Syntax:
%   val = fontGet(font, parm, [varargin])
%
% Description:
%    Get the value of the desired font parameter or property.
%
% Inputs:
%    font     - The font structure
%    param    - The parameter you wish to retrieve. Possibilites include:
%               type (Always 'font')
%               name
%               character
%               size
%               family
%               style
%               dpi
%               bitmap
%    varargin - (Optional) Additional data that may be required by a
%               parameter to retrieve an appropriate value.
%
% Outputs:
%    val      - The value of the specified parameter
%
% Optional key/value pairs:
%    **Needs to be filled out**
%
% Notes:
%    * [Note: JNM - Assign someone to fill out key/value pairs section.]
%

% History:
%    xx/xx/14  HJ/BW  PDCSoft Team, 2014.
%    02/27/18  jnm    Formatting

if notDefined('param'), error('Parameter must be defined.'); end

% Default is empty when the parameter is not yet defined.
val = [];

param = ieParamFormat(param);

switch param
    % Book keeping
    case 'type'
        val = font.type;
    case 'name'
        val = font.name;
    case 'character'
        val = font.character;
    case 'size'
        val = font.size;
    case 'family'
        val = font.family;
    case 'style'
        val = font.style;
    case 'dpi'
        val = font.dpi;
    case 'bitmap'
        val = font.bitmap;

    % Derived
    case 'paddedbitmap'
        % fontGet(font, 'padded bitmap', padval);
        % vcNewGraphWin;
        % imagesc(fontGet(font, 'padded bitmap'));
        % axis equal
        if ~isempty(varargin), padsize = varargin{1}; end
        if length(varargin) > 1, padval = varargin{2}; end

        if notDefined('padsize'), padsize = [7 7]; end
        if notDefined('padval'), padval = 1; end

        % RGB bitmap
        bitmap = fontGet(font, 'bitmap');

        % Pad and return
        newSize = size(bitmap); 
        newSize(1:2) = newSize(1:2) + 2 * padsize;
        val = zeros(newSize);        
        for ii = 1:size(bitmap, 3)
            val(:, :, ii) = padarray(bitmap(:, :, ii), padsize, padval);
        end

    otherwise
        disp(['Unknown parameter: ', param]);

 end

end