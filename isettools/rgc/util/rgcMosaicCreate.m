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
%   GLM     - Pillow et al. coupled generalized line
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
% See also: irCreate.m, rgcMosaic.m, rgcLinear.m, rgcLNP.m,
%           rgcGLM.m, t_rgc.m, t_rgcIntroduction.
%
% Copyright ISETBIO Team 2016
% 7/2016 JRG updated

%% Parse inputs

p = inputParser; 
p.addRequired('ir');

% Experiment ... thinking about input parsing more generally (JRG/BW)
mosaicTypes = {'onparasol','offparasol','onmidget','offmidget','smallbistratified','sbc'};
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
        % Chichilnisky & Kalmar, J. Neurosci (2002)
        obj = rgcLinear(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'lnp', 'rgclnp'}
        % Standard linear nonlinear poisson        
        % Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, J. Neurosci (2005);
        obj = rgcLNP(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case {'glm','rgcglm'}
        % Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli,
        % Nature (2008).
        obj = rgcGLM(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    case{'phys','rgcphys'}
        % Unit testing of the physiology
        obj = rgcPhys(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    otherwise
        error('Unknown inner retina class: %s\n',class(ir));
end
