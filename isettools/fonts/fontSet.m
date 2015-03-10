function val = fontSet(font, param, val, varargin)
%Get font parameters and derived properties
%
%     val = fontGet(font,parm,val,varargin)
%
%
% (HJ/BW)  PDCSoft Team, 2014.

if notDefined('param'), error('Parameter must be defined.'); end
if ~exist('val','var'), error('Parameter must be defined.'); end

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
        disp(['Unknown parameter: ',param]);
        
end

end
