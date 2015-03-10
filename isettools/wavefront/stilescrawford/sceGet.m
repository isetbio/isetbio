function val = sceGet(sceP,parm,varargin)
% val = sceGet(sceP,parm,varargin)
%
% Get function for Stiles-Crawford effect parameters
%
% Example:
%   sceGet(sceP,'xo')
%   sceGet(sceP,'wavelengths')
%   sceGet(sceP,'rho',550)
%
% See also sceCreate
%
% (c) Wavefront Toolbox Team, 2012

if notDefined('sceP'), error('sceP must be defined.'); end
if notDefined('parm'), error('Parameter must be defined.'); end

% Default is empty when the parameter is not yet defined.
parm = ieParamFormat(parm);

switch parm
    case 'xo'
        val = sceP.xo;
    case 'yo'
        val = sceP.yo;
    case {'wave','wavelengths'}
        val = sceP.wavelengths;
    case 'rho'
        if isempty(varargin)
            val = sceP.rho;
            return
        else
            wList = varargin{1};
            [l, loc] = ismember(wList,sceP.wavelengths); %#ok<ASGLU>
            val = sceP.rho(loc);
        end
    otherwise
        error('Unknown parameter %s\n',parm);
end

end
