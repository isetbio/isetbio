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

p.parse(obj,ir,varargin{:});

experimentID = p.Results.experimentID;
stimulusFit = p.Results.stimulusFit;
stimulusTest = p.Results.stimulusTest;
cellType =  p.Results.cellType;

obj.experimentID = experimentID;
obj.stimulusFit = stimulusFit;
obj.stimulusTest = stimulusTest;
obj.cellType = cellType;

%% Set defaults
obj.generatorFunction = @exp;
obj.numberTrials = 10;

% % Coupled experiment
% glmFitPath = '/Users/james/Documents/matlab/NSEM_data/';
% matFileNames = dir([glmFitPath '/ON*.mat']);

% switch ieParamFormat(stimulusFit)
%     case 'wn'        
%         
%         switch ieParamFormat(stimulusTest)
%             case 'wn'
%                 glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/';
%             case 'nsem'
%                 glmFitPath = '/Users/james/Documents/matlab/akheitman/WN_mapPRJ/Test_NSEM/';
%         end
%         
%     otherwise % case 'NSEM'
%         glmFitPath = '/Users/james/Documents/matlab/akheitman/NSEM_mapPRJ/';
% end



% RDT initialization
rdt = RdtClient('isetbio');
rdt.crp('resources/data/rgc');

switch ieParamFormat(cellType)
    case 'offparasol'
%         matFileNames = dir([glmFitPath experimentID '/OFF*.mat']);        
        data = rdt.readArtifact('mosaicGLM_WN_OFFParasol_2013_08_19_6', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;
        
        data2 = rdt.readArtifact('goodind_2013_08_19_6_OFFParasol', 'type', 'mat');
        goodind = data2.goodind;
    otherwise % case 'onparasol'
%         matFileNames = dir([glmFitPath experimentID '/ON*.mat']);
        data = rdt.readArtifact('mosaicGLM_WN_ONParasol_2013_08_19_6', 'type', 'mat');
        mosaicGLM = data.mosaicGLM;        
        
        data2 = rdt.readArtifact('goodind_2013_08_19_6_ONParasol', 'type', 'mat');
        goodind = data2.goodind;
end

goodind = 1:length(mosaicGLM);
                        
% % % % % % 
% Loop through mat files and load parameters
% for matFileInd = 1:length(mosaicGLM)
for matFileInd = 1:31%length(goodind)
%     cell = matFileNames(matFileInd).name(1:end-4);
%     obj.cellID{matFileInd,1} = cell;
%     load([glmFitPath experimentID '/' cell '.mat']);

    obj.cellID{matFileInd,1} = mosaicGLM{matFileInd}.cell_savename;
%     obj.cellID{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.cell_savename;
    
    obj.postSpikeFilter{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.PostSpike.Filter;
    if isfield(mosaicGLM{goodind(matFileInd)}.linearfilters,'Coupling')

        obj.couplingFilter{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Coupling.Filter;
    end
    
    obj.tonicDrive{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.TonicDrive.Filter;
    
    obj.sRFcenter{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.space_rk1;
    obj.sRFsurround{matFileInd,1} = 0*mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.space_rk1;
    obj.tCenter{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.time_rk1;
    obj.tSurround{matFileInd,1} = 0*mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.time_rk1;
    
    if isfield(mosaicGLM{goodind(matFileInd)}.linearfilters,'Coupling')

        couplingMatrixTemp{matFileInd,1} = mosaicGLM{goodind(matFileInd)}.cellinfo.pairs;
    end
    
  obj.cellLocation{matFileInd,1} = [mosaicGLM{goodind(matFileInd)}.cellinfo.slave_centercoord.x_coord mosaicGLM{goodind(matFileInd)}.cellinfo.slave_centercoord.y_coord];
    
end

obj.rfDiameter = size(mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.Filter,1);