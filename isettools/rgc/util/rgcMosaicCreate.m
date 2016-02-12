function ir = rgcMosaicCreate(ir, varargin)
%% Add a specific type of RGC mosaic to an inner retina (IR) object
%
%  ir = rgcMosaicCreate(ir, mosaicType)
%
% We can add these types of mosaics to the inner retina object:
%
%       ON Parasol
%       OFF Parasol
%       ON Midget
%       OFF Midget
%       Small Bistratified
%
% Each one of these types will be implemented for the various computational
% models (e.g., GLM, LNP and so forth).  The model is defined before the
% creation of the IR object, and the class of the IR object is dependent on
% the model.
% 
% This function can alternatively be called using the object.method()
% syntax.
% 
% Examples:
%
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

