function ir = rgcMosaicCreate(ir, varargin)
% Add a specific type of RGC mosaic to a specific type of computation 
%
%  ir = rgcMosaicCreate(ir, mosaicType)
%
% The computation (e.g, GLM,Linear,LNP) is specified by the class of the
% inner retina (ir)
%
% The implemented mosaic types are
%
% {'ON Parasol', 'OFF Parasol', 'ON Midget', 'OFF Midget', 'Small
% Bistratified' }
%
% Examples:
%       innerRetina = rgcMosaicCreate(innerRetina,'mosaicType','on parasol');
%       innerRetina.mosaicCreate('mosaicType','on midget');
% 
% See also: rgcMosaic.m, rgcMosaicLinear.m, rgcMosaicLNP.m, rgcMosaicGLM.m,
%               t_rgc.m, t_rgcIntroduction.
%
% Copyright ISETBIO Team 2016

%% Parse inputs

p = inputParser; 
p.addRequired('ir');

mosaicTypes = {'on parasol','off parasol','on midget','off midget','small bistratified','sbc'};
p.addParameter('mosaicType','on parasol',@(x) any(validatestring(x,mosaicTypes)));

modelTypes = {'linear','lnp','glm','phys','subunit','pool'};
p.addParameter('model','linear',@(x) any(validatestring(x,modelTypes)));

p.parse(ir,varargin{:});

%% Specify the ganglion cell mosaic type
mosaicType = p.Results.mosaicType;
model = p.Results.model;
%% Switch on the computational model

% There is a separate class for each ir computational model.  These are
% rgcglm, rgclinear ...
switch ieParamFormat(model)
    case {'linear','rgclinear'}
        obj = rgcMosaicLinear(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'pool', 'rgcpool'}
        obj = rgcMosaicPool(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'lnp', 'rgclnp'}
        obj = rgcMosaicLNP(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'glm','rgcglm'}
        obj = rgcMosaicGLM(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'subunit','rgcsubunit'}
        obj = rgcMosaicSubunit(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case{'phys','rgcphys'}
        obj = rgcMosaicPhys(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    otherwise
        error('Unknown inner retina class: %s\n',class(ir));
end

