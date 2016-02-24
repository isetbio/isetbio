function ir = rgcMosaicCreate(ir, varargin)
% Add a type of RGC mosaic with a specific computational model
%
% The rgc mosaics are stored as a cell array within the inner retina (ir)
% class.  The RGC mosaics are the main computational engine for producing
% retinal spike outputs.
%
% The implemented computational models are 'linear', 'LNP', or 'GLM'
%
% The implemented mosaic types are
%   'ON Parasol', 
%   'OFF Parasol', 
%   'ON Midget', 
%   'OFF Midget', 
%   'Small Bistratified' 
%
% Examples:
%  This function can be called as:
%
%   ir = rgcMosaicCreate(ir, 'model', ['linear','GLM',etc.], 'mosaicType', ['on parasol', 'sbc', etc.])
%
%  Often, we call it as a method of the inner retina class. In that case,
%  the call looks like:
%
%   ir = irCreate(osCreate('identity'));
%   ir.mosaicCreate('model','linear','mosaicType','on parasol');
%   ir.mosaicCreate('model','GLM','mosaicType','on midget');
%
% See also: irCreate, rgcMosaic.m, rgcMosaicLinear.m, rgcMosaicLNP.m, rgcMosaicGLM.m,
%               t_rgc.m, t_rgcIntroduction.
%
% Copyright ISETBIO Team 2016

%% Parse inputs

p = inputParser; 
p.addRequired('ir');

% Experiment ... thinking about input parsing more generally (JRG/BW)
mosaicTypes = {'on parasol','off parasol','on midget','off midget','small bistratified','sbc'};
% p.addParameter('mosaicType','on parasol',@(x) any(validatestring(x,mosaicTypes)));
p.addParameter('type','on parasol',@(x) any(validatestring(x,mosaicTypes)));
modelTypes = {'linear','lnp','glm','phys','subunit','pool'};
p.addParameter('model','linear',@(x) any(validatestring(x,modelTypes)));

p.parse(ir,varargin{:});

%% Specify the ganglion cell mosaic type
mosaicType = p.Results.type;
model = p.Results.model;
%% Switch on the computational model

% There is a separate mosaic class for each ir computational model.  
% These are rgcMosaicLinear, rgcMosaicLNP, rgcMosaicGLM,...
switch ieParamFormat(model)
    case {'linear','rgclinear'}
        obj = rgcLinear(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'pool', 'rgcpool'}
        obj = rgcPool(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'lnp', 'rgclnp'}
        obj = rgcLNP(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'glm','rgcglm'}
        obj = rgcGLM(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'subunit','rgcsubunit'}
        obj = rgcSubunit(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case{'phys','rgcphys'}
        obj = rgcPhys(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    otherwise
        error('Unknown inner retina class: %s\n',class(ir));
end

