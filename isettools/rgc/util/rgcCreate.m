function obj = rgcCreate(varargin)
% rgcCreate: generate an @rgcLinear, @rgcLNP or @rgcGLM object.
% 
% Inputs: none.
% Outputs: the rgc object.
% 
% See the initialize method for the @rgcLinear, @rgcLNP or @rgcGLM 
% subclasses for more details of the specific implementations.
%
% 9/2015 JRG

if nargin ~= 7
        warning(sprintf(['\nrgcCreate: create an isetbio @rgc object.\n'...
            'rgcCreate(type, sensor, outersegment, eyeSide, patchRadius, patchAngle)\n'...
            'rgc1 = rgcCreate(''Linear'', scene, sensor, identityOS, ''right'', 3.0, 180);']));
        return;
end

rgcType = 'linear';    % Default values
if ~isempty(varargin),   rgcType = ieParamFormat(varargin{1}); end

%% Create the proper object
switch rgcType
    case {'linear','rgclinear'}
        obj = rgcLinear(varargin{2:7});
    case {'lnp','rgclnp'}
        obj = rgcLNP(varargin{2:7});
    case {'glm','rgcglm'}
        obj = rgcGLM(varargin{2:7});
    case {'subunit','rgcsubunit'}
        obj = rgcSubunit(varargin{2:7});
    otherwise
        obj = rgcLinear(varargin{2:7});
end
