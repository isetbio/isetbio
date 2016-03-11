function glmprs = setGLMprs(mosaic)
% Converts the mosaic object into the appropriate format for the Pillow
% spike-generating code.
% 
% The properties of 'glmprs' from "testscript_GLM_coupled.m" found in
% \pillow_code_GLM_v1_Feb2010, available at 
% http://pillowlab.princeton.edu/code_GLM.html.
% 
% Target output structure:
% glmprs = 
% 
%            k: [20x1x2 double]
%           dc: [3 3]
%           ih: [598x2x2 double]
%          iht: [598x1 double]
%           dt: 0.0100
%        nlfun: @exp
%     ihbasprs: [1x1 struct]
% 
% 3/2016 JRG (c) isetbio

%% Get number of cells
nCells = size(mosaic.cellLocation);

%% Set linear temporal filter
% Because of the construction of simGLM and simGLMcpl, we compute the
% linear spatial and temporal responses before passing to simGLM. In order
% to not do any temporal filtering internal to simGLM, we set the temporal
% responses to impulses. The convolution of each neuron's linear response
% with an impulse results in the original linear response.

glmprs.k = zeros(length(mosaic.tCenter{1}),1,2);
% glmprs.k(:,1,1) = mosaic.tCenter{1};
% glmprs.k(:,1,2) = mosaic.tSurround{1};
glmprs.k = zeros(20,nCells(1)*nCells(2),nCells(1)*nCells(2));

% Set impulse
cellCtr = 0;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
        glmprs.k(20,cellCtr,cellCtr) = 1;
    end
end

%% Set DC value
glmprs.dc = zeros(size(mosaic.sRFcenter));

%% Set coupling filters
ihcpl = mosaicGet(mosaic, 'couplingFilter');

hlen = length(ihcpl{1,1});
nCells = size(mosaic.cellLocation);
nCellsTotal = nCells(1)*nCells(2);
spikeTimes = cell(nCells);
cellCtr = 0;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
        ih(cellCtr,:,:) = reshape(mosaic.couplingFilter{ind2sub([nCells(1),nCells(2)],cellCtr)},nCellsTotal,hlen);
%         ih(cellCtr,:,:) = reshape(mosaic.couplingFilter{ind2sub([nCells(1),nCells(2)],cellCtr)},nCellsTotal,hlen);

    end
end

% [sz1,sz2,sz3]=size(ih);
% ih = reshape(ih,[sz1 sz
    
ih = permute(ih,[3 2 1]); % flip 2nd & 3rd dimensions

glmprs.ih = ih;

%% Set time samples
glmprs.iht = mosaic.dt*(1:length(ih(:,1,1)));

%% Set interpolation
glmprs.dt = mosaic.dt;

%% Set nonlinearity
glmprs.nlfun = mosaic.generatorFunction;



