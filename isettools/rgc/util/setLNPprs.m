function glmprs = setGLMprs(mosaic)
% Convert a mosaic object to a pillow spike-generating code format.
%
% Syntax:
%   glmprs = setGLMprs(mosaic)
%
% Description:
%    Converts the mosaic object into the appropriate format for the Pillow
%    spike-generating code.
% 
%    Target output structure:
%       ggsim = 
%            type: 'glm'
%               k: [20x1 double]
%           nlfun: @exp
%              dc: 3
%              ih: [598x1 double]
%             iht: [598x1 double]
%              dt: 0.0100
%         kbasprs: []
%        ihbasprs: [1x1 struct]
%
% Inputs:
%    mosaic - Object. An rgc mosaic object.
%
% Outputs:
%    glmprs - Struct. A structure in a format that is applicable for the
%             pillow spike-generating code. The format of the structure can
%             be found above.
%
% Optional key/value pairs:
%    None.
%
% References:
%    * The properties of 'glmprs' from "testscript_GLM_coupled.m" found in
%      \pillow_code_GLM_v1_Feb2010, available at 
%      http://pillowlab.princeton.edu/code_GLM.html.
%

% History:
%    03/XX/16  JRG  (c) isetbio
%    05/31/19  JNM  Documentation pass

%% Get number of cells
nCells = size(mosaic.cellLocation);

%% Set linear temporal filter
% Because of the construction of simGLM and simGLMcpl, we compute the
% linear spatial and temporal responses before passing to simGLM. In order
% to not do any temporal filtering internal to simGLM, we set the temporal
% responses to impulses. The convolution of each neuron's linear response
% with an impulse results in the original linear response.

glmprs.k = zeros(length(mosaic.tCenter{1}), 1, 2);
% glmprs.k(:, 1, 1) = mosaic.tCenter{1};
% glmprs.k(:, 1, 2) = mosaic.tSurround{1};
glmprs.k = zeros(20, nCells(1) * nCells(2), nCells(1) * nCells(2));

% Set impulse
cellCtr = 0;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr + 1;
        glmprs.k(20, cellCtr, cellCtr) = 1;
    end
end

%% Set DC value
glmprs.dc = zeros(size(mosaic.sRFcenter));

%% Set coupling filters
% For only PSF, set all coupling filters to zero
postSpikeFilter = mosaicGet(mosaic, 'postSpikeFilter');

hlen = length(postSpikeFilter);
nCells = size(mosaic.cellLocation);
nCellsTotal = nCells(1) * nCells(2);
spikeTimes = cell(nCells);
cellCtr = 0;
ih = zeros(nCellsTotal, nCellsTotal, 598);
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr + 1;
        if xcell == ycell
            ih(cellCtr, cellCtr, :) = zeros(1, 598);
            %reshape(postSpikeFilter, nCellsTotal, hlen);
        else
            ih(cellCtr, :, :) = zeros(nCellsTotal, 598);
        end
    end
end

% [sz1, sz2, sz3]=size(ih);
% ih = reshape(ih, [sz1 sz

ih = permute(ih, [3 2 1]); % flip 2nd & 3rd dimensions

glmprs.ih = ih;
%% Set time samples
glmprs.iht = mosaic.dt * (1:length(glmprs.ih));

%% Set interpolation
glmprs.dt = mosaic.dt;

%% Set nonlinearity
glmprs.nlfun = mosaic.generatorFunction;
