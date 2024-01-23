function val = conerectGet(cm,param,varargin)
% Get method for rectangular cone mosaics
%
% Synopsis
%   val = conrectGet(cm,param,varargin)
%
% Brief description
%   Many useful parameters (dependent variables) derived from the
%   coneRectMosaic class
%
%   Not sure we should do it this way or using get.variable name.
%

if notDefined('cm'), error('coneRectMosaic required.'); end
if notDefined('param'), error('Parameter required.'); end
val = [];

switch param
    case 'absorptions'
        % cm.get('absorptions','conetype')
        if isempty(varargin)
            val = cm.absorptions;
        else
            switch ieParamFormat(varargin{1})
                case {'lcones','l'}
                    val = cm.absorptions(cm.pattern == 2);
                case {'mcones','m'}
                    val = cm.absorptions(cm.pattern == 3);
                case {'scones','s'}
                    val = cm.absorptions(cm.pattern == 4);
            end
        end

    otherwise
        error('Unknown parameter %s\n',param);
end

end

