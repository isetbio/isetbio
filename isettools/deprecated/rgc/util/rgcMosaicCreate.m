function rgcM = rgcMosaicCreate(rgcL, bipolarM, cellType, varargin)
% Add an RGC mosaic with a specific computational model and cell type
%
% Syntax:
%   rgcM = rgcMosaicCreate(rgcL, bipolarM, cellType, [varargin])
%
% Description:
%    The rgc mosaics are stored as a cell array attached to the an rgcLayer
%    class.  The rgc mosaics are the computational engine for producing
%    retinal spike outputs.
%
%    The rgc mosaic types are
%        'ON Parasol', 
%        'OFF Parasol', 
%        'ON Midget', 
%        'OFF Midget', 
%        'ON sbc' 
%
%    The  models are
%        GLM     - Pillow et al. coupled generalized line
%        LNP     - Linear, nonlinear, poisson (EJ 2002 reference)
%        Phys    - Fitting the physiology data from EJ
%
%    Often, we call it as a method of the rgcLayer class. In that case, 
%    the call looks like: 
%         rgcL.mosaicCreate('model', val, 'type', val);
%
% Inputs:
%    rgcL     - Object. An RGC Layer object.
%    bipolarM - Object. A bipolar mosaic object.
%    cellType - String. A string describing the rgc mosaic type to create.
%               The options are: 'onParasol', 'offParasol', 'onMidget',
%               'offMidget', and 'onSBC'
%
% Outputs:
%    rgcM     - Object. An RGC Mosaic object.
%
% Optional key/value pairs:
%    model    - String. A string describing the model type. Default 'lnp'.
%               The options are: 'linear', 'lnp', 'glm', 'phys', 'subunit',
%               and 'pool'.
%
% See Also:
%   s_LayersTest.m
%

% History:
%    XX/XX/16       Copyright ISETBIO Team 2016
%    05/31/19  JNM  Documentation pass

%% Parse inputs
p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('rgcL', @(x)(isequal(class(x), 'rgcLayer')));
p.addRequired('bipolarM', @(x)(isequal(class(x), 'bipolarMosaic')));

% Experiment ... thinking about input parsing more generally (JRG/BW)
cellTypes = {'onparasol', 'offparasol', 'onmidget', 'offmidget', 'onsbc'};
p.addRequired('cellType', ...
    @(x) any(validatestring(ieParamFormat(x), cellTypes)));

modelTypes = {'linear', 'lnp', 'glm', 'phys', 'subunit', 'pool'};
p.addParameter('model', 'lnp', ...
    @(x) any(validatestring(ieParamFormat(x), modelTypes)));

p.parse(rgcL, bipolarM, cellType, varargin{:});
model = p.Results.model;

%% Switch on the computational model
% There is a separate mosaic class for each ir computational model.  
% These are rgcMosaicLinear, rgcMosaicLNP, rgcMosaicGLM, ...
switch ieParamFormat(model)
    case {'lnp'}
        % Standard linear nonlinear poisson
        % Pillow, Paninski, Uzzell, Simoncelli & Chichilnisky, ...
        % J. Neurosci (2005);
        rgcM = rgcLNP(rgcL, bipolarM, cellType, p.Unmatched);
    case {'glm'}
        % Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky & ...
        % Simoncelli, Nature (2008).
        rgcM = rgcGLM(rgcL, bipolarM, cellType, p.Unmatched);
    case{'phys'}
        % Unit testing of the physiology
        % Requires the isetbio repository EJLExperimentalRGC
        rgcM = rgcPhys(rgcL, bipolarM, cellType);
    otherwise
        error('Unknown inner retina class: %s\n', class(rgcL));
end

%% We want to selecct a particular mosaic with nMosaic here
% As written, add this mosaic to the end of the list.
rgcL.set('mosaic', rgcM);

end
