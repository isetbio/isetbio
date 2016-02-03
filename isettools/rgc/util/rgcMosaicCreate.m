function rgc = rgcMosaicCreate(rgc, varargin)
% rgcCreate: generate an rgcMosaicLinear, rgcMosaicLNP or rgcMosaicGLM object.
%
%  obj = rgcCreate(rgc, mosaicType)
% 


narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% addRequired(p,'rgc',@isa();
addParameter(p,'mosaicType','onParasol',@ischar);

p.parse(varargin{:});

% Set key-value pairs
switch ieParamFormat(class(rgc))
    case {'linear','rgclinear'}
        obj = rgcMosaicLinear(p.rgc);
        rgcSet(rgc, 'mosaic', rgc);
    case {'pool','rgcpool'}
        obj = rgcMosaicPool(rgc);
        rgcSet(rgc, 'mosaic', obj);
    case {'lnp','rgclnp'}
        obj = rgcMosaicLNP(rgc);
        rgcSet(rgc, 'mosaic', obj);
    case {'glm','rgcglm'}
        obj = rgcMosaicGLM(rgc);
        rgcSet(rgc, 'mosaic', obj);
    case {'subunit','rgcsubunit'}
        obj = rgcMosaicSubunit(rgc);
        rgcSet(rgc, 'mosaic', obj);
    case{'phys','rgcphys'}
        obj = rgcMosaicPhys(rgc);
        rgcSet(rgc, 'mosaic', obj);
    otherwise
        obj = rgcMosaicLinear(rgc);
        rgcSet(rgc, 'mosaic', obj);        
end

