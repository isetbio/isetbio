function obj = initialize(obj, ir, varargin)
% Initializes the rgcPhys object by loading a mosaic of GLM fits from an
% experiment in the Chichilnisky lab.
% 
% This function is only called by rgcPhys, which itself is only called by irPhys.
% 
%             rgcPhys = rgcPhys.initialize(rgc, varargin{:});   
% Inputs: 
%       rgc: an isetbio rgcPhys object
% 
% Outputs: the mosaic object, where each cell has a location, linear spatial
% and temporal receptive fields, a DC offest, a generator function, a
% post-spike filter, coupling filters if necessary, and empty fields for
% the linear, voltage and spiking responses.
% 
% 
% See also rgcPhys, irPhys.
% 
% (c) isetbio
% 09/2015 JRG%


%% Parse inputs
p = inputParser;
p.addRequired('obj');
p.addRequired('ir');
addParameter(p,'experimentID','2013-08-19-6',@ischar);
addParameter(p,'stimulusFit','WN',@ischar);
addParameter(p,'stimulusTest','NSEM',@ischar);
addParameter(p,'cellType','OnParasol',@ischar);


addParameter(p,'name','inner retina',@ischar);
addParameter(p,'species','unknown',@ischar);
% addParameter(p,'outersegment','os',@ischar);
addParameter(p,'eyeSide',    'left', @ischar);
addParameter(p,'eyeRadius',   4,     @isnumeric);
addParameter(p,'eyeAngle',    0,     @isnumeric);  % X-axis is 0, positive Y is 90
addParameter(p,'fov',         8,     @isnumeric);
addParameter(p,'cellIndices',   0,     @isnumeric);
addParameter(p,'inputSize',80,           @isnumeric);
addParameter(p,'inputScale',1,           @isnumeric);
addParameter(p,'averageMosaic', 0,      @isnumeric);

p.parse(obj,ir,varargin{:});

experimentID = p.Results.experimentID;
stimulusFit  = p.Results.stimulusFit;
stimulusTest = p.Results.stimulusTest;
cellType     =  p.Results.cellType;
cellIndices  = p.Results.cellIndices;
inputSize    = p.Results.inputSize;
inputScale   = p.Results.inputScale;

fov = p.Results.fov;
ecc = p.Results.eyeRadius;

averageFlag = p.Results.averageMosaic;

obj.experimentID = experimentID;
obj.stimulusFit  = stimulusFit;
obj.stimulusTest = stimulusTest;
obj.cellType     = cellType;

%% Set defaults
% obj.generatorFunction = @exp;
obj.numberTrials = 10;

[mosaicGLM, goodind] = glmLoad(obj, 'cellType', cellType);

switch(averageFlag)
    
    % Load in a mosaic from a physiology experiment
    case 0        
        obj = mosaicLoadExperimental(obj, mosaicGLM, cellType, cellIndices, goodind);
        
    % Generate a mosaic that fully tiles the visual field 
    % with the average properties of an experimental mosaic,
    % shifted to a specified eccentricity.
    case 1        
        obj = mosaicLoadAverage(obj, mosaicGLM, cellType, 'cellIndices',cellIndices, 'goodind',goodind, 'params', p.Results);
end