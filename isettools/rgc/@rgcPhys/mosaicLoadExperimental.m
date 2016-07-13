function obj = mosaicLoadExperimental(obj, mosaicGLM, cellType, varargin)
% Load the parameters for an RGC mosaic measured in an experiment by the
% Chichilnisky Lab.
% 
% JRG (c) 2016 isetbio team

%% Parse inputs
p = inputParser;
p.addRequired('obj');
p.addRequired('mosaicGLM');
p.addRequired('cellType');
addParameter(p,'cellIndices',   4,     @isnumeric);
addParameter(p,'goodind',    0,     @isnumeric);  

p.parse(obj,mosaicGLM,cellType,varargin{:});

cellType = p.Results.cellType;
cellIndices = p.Results.cellIndices;
goodind  = p.Results.goodind;

%% Select which cells to load
if cellIndices ~= 0
    cellIndicesEval = cellIndices;
else
    cellIndicesEval = [1:length(goodind)];
end
matFileCtr = 0;

%% Loop through mosaic GLM and load parameters

for matFileInd = cellIndicesEval
    
    matFileCtr = matFileCtr+1;
    obj.cellID{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.cell_savename;
  
    % Post spike filter
    if isfield(mosaicGLM{goodind(matFileInd)}.linearfilters,'PostSpike')
        obj.postSpikeFilter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.PostSpike.Filter;
    else
        obj.postSpikeFilter{matFileCtr,1} = 0;
    end
    
    % Coupling filter
    if isfield(mosaicGLM{goodind(matFileInd)}.linearfilters,'Coupling')
        obj.couplingFilter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Coupling.Filter;
         couplingMatrixTemp{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.cellinfo.pairs;
    end
        
    % Tonic drive 
    switch ieParamFormat(cellType)
        case {'onparasolrpe','offparasolrpe','onmidgetrpe','offmidgetrpe','onsbcrpe','sbcrpe'}
            obj.tonicDrive{matFileCtr,1} = 0;
        otherwise
            obj.tonicDrive{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.TonicDrive.Filter;
    end
    
    % Generator function
    switch ieParamFormat(cellType)
        case {'onparasolrpe','offparasolrpe','onmidgetrpe','offmidgetrpe','onsbcrpe','sbcrpe'}
            obj.generatorFunction{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.model;
        otherwise
            obj.generatorFunction{matFileCtr,1} = @exp;
    end
    
    % Spatial receptive field
    obj.sRFcenter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.space_rk1;
    obj.sRFsurround{matFileCtr,1} = 0*mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.space_rk1;
    obj.tCenter{matFileCtr,1} = mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.time_rk1;
    obj.tSurround{matFileCtr,1} = 0*mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.time_rk1;
    
    % Spatial receptive field location
    obj.cellLocation{matFileCtr,1} = [mosaicGLM{goodind(matFileInd)}.cellinfo.slave_centercoord.x_coord mosaicGLM{goodind(matFileInd)}.cellinfo.slave_centercoord.y_coord];
    
end

% RF size in stimulus pixels, same for every cell
obj.rfDiameter = size(mosaicGLM{goodind(matFileInd)}.linearfilters.Stimulus.Filter,1);
