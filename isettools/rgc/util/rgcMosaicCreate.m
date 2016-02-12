function ir = rgcMosaicCreate(ir, varargin)
% Add a specific type of RGC mosaic to an inner retina (ir) object
%
%  ir = rgcCreate(ir, mosaicType)
%
% We can add these types of mosaics to the inner retina object
%
%   on/off parasol
%   on/off midget
%   small bistratified
%
% Each one of these types will be implemented for the various computational
% models (e.g., GLM, LNP and so forth).  The model is a parameter in the ir
% object itself.
% 
% Example:
%
% See also
%
% Copyright ISETBIO Team 2016

%% Parse inputs

p = ieInputParser; 
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

