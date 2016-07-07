function ir = rgcMosaicCreate(ir, varargin)
% Add a type of RGC mosaic with a specific computational model
%
% The rgc mosaics are stored as a cell array within the inner retina (ir)
% class.  The RGC mosaics are the main computational engine for producing
% retinal spike outputs.
%
% The implemented mosaic types are
%
%   'ON Parasol', 
%   'OFF Parasol', 
%   'ON Midget', 
%   'OFF Midget', 
%   'Small Bistratified' 
%
% The implemented models are
%
%   Linear  - Straight linear convolution, no spikes
%   GLM     - Pillow et al. coupled generalized linear model
%   Pool    - Built by Winawer for psychophysics
%   Subunit - Related to Markus Meister modeling
%   LNP     - Linear, nonlinear, poisson (EJ 2002 reference)
%   Phys    - Fitting the physiology data from EJ
%
% Examples:
%
%   ir = rgcMosaicCreate(ir, 'model', ['linear','GLM',etc.], 'mosaicType', ['on parasol', 'sbc', etc.])
%
% Often, we call it as a method of the inner retina class. In that case,
% the call looks like: 
%
%   ir = irCreate(osCreate('identity'));
%   ir.mosaicCreate('model','linear','type','on parasol');
%   ir.mosaicCreate('model','GLM','type','on midget');
%
% See also: irCreate, rgcMosaic.m, rgcMosaicLinear.m, rgcMosaicLNP.m,
%           rgcMosaicGLM.m, t_rgc.m, t_rgcIntroduction.
%
% Copyright ISETBIO Team 2016

%% Parse inputs

p = inputParser; 
p.addRequired('ir');

% Experiment ... thinking about input parsing more generally (JRG/BW)
mosaicTypes = {'onparasol','offparasol','onmidget','offmidget','smallbistratified','sbc'};
% p.addParameter('mosaicType','on parasol',@(x) any(validatestring(x,mosaicTypes)));
p.addParameter('type','on parasol',@(x) any(validatestring(ieParamFormat(x),mosaicTypes)));
modelTypes = {'linear','lnp','glm','phys','subunit','pool'};
p.addParameter('model','linear',@(x) any(validatestring(x,modelTypes)));

p.parse(ir,varargin{:});

%% Specify the ganglion cell mosaic type
mosaicType = p.Results.type;
model      = p.Results.model;
%% Switch on the computational model

% There is a separate mosaic class for each ir computational model.  
% These are rgcMosaicLinear, rgcMosaicLNP, rgcMosaicGLM,...
switch ieParamFormat(model)
    case {'linear','rgclinear'}
        % Straight linear convolution, no spikes
        obj = rgcLinear(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'pool', 'rgcpool'}
        % Built by Winawer for psychophysics
        obj = rgcPool(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'lnp', 'rgclnp'}
        % Standard linear nonlinear poisson
        % EJ 2002 reference
        obj = rgcLNP(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'glm','rgcglm'}
        % Pillow et al. 2008
        obj = rgcGLM(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'subunit','rgcsubunit'}
        % Related to Markus Meister modeling
        obj = rgcSubunit(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case{'phys','rgcphys'}
        % Unit testing of the physiology
        obj = rgcPhys(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    otherwise
        error('Unknown inner retina class: %s\n',class(ir));
end

