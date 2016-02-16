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
p.parse(ir,varargin{:});

%% Specify the ganglion cell mosaic type
mosaicType = p.Results.mosaicType;

%% Switch on the computational model

% There is a separate class for each ir computational model.  These are
% rgcglm, rgclinear ...
switch ieParamFormat(class(ir))
    case {'rgclinear'}
        obj = rgcMosaicLinear(ir, mosaicType);
        rgcSet(ir, 'mosaic', obj);
    case {'rgcpool'}
        obj = rgcMosaicPool(ir, mosaicType);
        rgcSet(ir, 'mosaic', obj);
    case {'rgclnp'}
        obj = rgcMosaicLNP(ir, mosaicType);
        rgcSet(ir, 'mosaic', obj);
    case {'rgcglm'}
        obj = rgcMosaicGLM(ir, mosaicType);
        rgcSet(ir, 'mosaic', obj);
    case {'rgcsubunit'}
        obj = rgcMosaicSubunit(ir, mosaicType);
        rgcSet(ir, 'mosaic', obj);
    case{'rgcphys'}
        obj = rgcMosaicPhys(ir, mosaicType);
        rgcSet(ir, 'mosaic', obj);
    otherwise
        error('Unknown inner retina class: %s\n',class(ir));
end

