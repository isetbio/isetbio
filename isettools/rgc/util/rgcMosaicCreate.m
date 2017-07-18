function rgcM = rgcMosaicCreate(rgcL, varargin)
% Add an RGC mosaic with a specific computational model and cell type
% 
%      rgcM = rgcMosaicCreate(rgcL, 'model',val, 'mosaicType',val)
% or
%
% The rgc mosaics are stored as a cell array attached to the an rgcLayer
% class.  The rgc mosaics are the computational engine for producing
% retinal spike outputs.
%
% The rgc mosaic types are
%     'ON Parasol', 
%     'OFF Parasol', 
%     'ON Midget', 
%     'OFF Midget', 
%     'ON sbc' 
%
% The  models are
%   GLM     - Pillow et al. coupled generalized line
%   LNP     - Linear, nonlinear, poisson (EJ 2002 reference)
%   Phys    - Fitting the physiology data from EJ
%
% Often, we call it as a method of the rgcLayer class. In that case,
% the call looks like: 
%
%      rgcL.mosaicCreate('model', val, 'type', val);
%
% See also: s_LayersTest.m
%
% Copyright ISETBIO Team 2016

%% Parse inputs

p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('ir');

% Experiment ... thinking about input parsing more generally (JRG/BW)
mosaicTypes = {'onparasol','offparasol','onmidget','offmidget','smallbistratified','onsbc'};
p.addParameter('type','on parasol',@(x) any(validatestring(ieParamFormat(x),mosaicTypes)));

modelTypes = {'linear','lnp','glm','phys','subunit','pool'};
p.addParameter('model','lnp',@(x) any(validatestring(ieParamFormat(x),modelTypes)));

p.parse(rgcL,varargin{:});

mosaicType = p.Results.type;
model      = p.Results.model;
%% Switch on the computational model

% There is a separate mosaic class for each ir computational model.  
% These are rgcMosaicLinear, rgcMosaicLNP, rgcMosaicGLM,...
switch ieParamFormat(model)
    case {'lnp', 'rgclnp'}
        % Standard linear nonlinear poisson        
        % Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, J. Neurosci (2005);
        rgcM = rgcLNP(rgcL, mosaicType,p.Unmatched);
        irSet(rgcL, 'mosaic', rgcM);
    case {'glm','rgcglm'}
        % Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & Simoncelli,
        % Nature (2008).
        rgcM = rgcGLM(rgcL, mosaicType,p.Unmatched);
        irSet(rgcL, 'mosaic', rgcM);
    case{'phys','rgcphys'}
        % Unit testing of the physiology
        % Requires the isetbio repository EJLExperimentalRGC
        rgcM = rgcPhys(rgcL, mosaicType);
        irSet(rgcL, 'mosaic', rgcM);
    otherwise
        error('Unknown inner retina class: %s\n',class(rgcL));
end

end
