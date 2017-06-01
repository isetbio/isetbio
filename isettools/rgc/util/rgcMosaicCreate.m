function ir = rgcMosaicCreate(ir, varargin)
% Add an RGC mosaic with a specific computational model and cell type
% 
%      ir = rgcMosaicCreate(ir, 'model',val, 'mosaicType',val)
% or
%      ir.mosaicCreate('model', val, 'type', val);
%
% The rgc mosaics are stored as a cell array attached to the inner retina
% (ir) class.  The RGC mosaics are the main computational engine for
% producing retinal spike outputs.
%
% The  mosaic types are
%     'ON Parasol', 
%     'OFF Parasol', 
%     'ON Midget', 
%     'OFF Midget', 
%     'Small Bistratified' 
%
% The  models are
%   GLM     - Pillow et al. coupled generalized line
%   LNP     - Linear, nonlinear, poisson (EJ 2002 reference)
%   Phys    - Fitting the physiology data from EJ
%
% Often, we call it as a method of the inner retina class. In that case,
% the call looks like: 
%
%     ir = irCreate(bipolar(coneMosaic));
%     ir.mosaicCreate('model','lnp','type','on parasol');
%     ir.mosaicCreate('model','GLM','type','on midget');
%
% See also:   t_rgc<>.m,     v_rgc<>, contain many examples.  These may
% rely on the RemoteData Toolbox downloads of pre-computed stimuli.
%
% Copyright ISETBIO Team 2016

%% Parse inputs

p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('ir');

% Experiment ... thinking about input parsing more generally (JRG/BW)
mosaicTypes = {'onparasol','offparasol','onmidget','offmidget','smallbistratified','sbc'};
p.addParameter('type','on parasol',@(x) any(validatestring(ieParamFormat(x),mosaicTypes)));

modelTypes = {'linear','lnp','glm','phys','subunit','pool'};
p.addParameter('model','lnp',@(x) any(validatestring(x,modelTypes)));

p.parse(ir,varargin{:});

mosaicType = p.Results.type;
model             = p.Results.model;
%% Switch on the computational model

% There is a separate mosaic class for each ir computational model.  
% These are rgcMosaicLinear, rgcMosaicLNP, rgcMosaicGLM,...
switch ieParamFormat(model)
    case {'lnp', 'rgclnp'}
        % Standard linear nonlinear poisson        
        % Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, J. Neurosci (2005);
        obj = rgcLNP(ir, mosaicType,p.Unmatched);
        irSet(ir, 'mosaic', obj);
    case {'glm','rgcglm'}
        % Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli,
        % Nature (2008).
        obj = rgcGLM(ir, mosaicType,p.Unmatched);
        irSet(ir, 'mosaic', obj);
    case{'phys','rgcphys'}
        % Unit testing of the physiology
        % Requires the isetbio repository EJLExperimentalRGC
        obj = rgcPhys(ir, mosaicType);
        irSet(ir, 'mosaic', obj);
    otherwise
        error('Unknown inner retina class: %s\n',class(ir));
end

end
