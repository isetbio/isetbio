function spikePrs = setSpikePrs(mosaic)
% Use the pillow spike-generating code to convert the mosaic format.
%
% Syntax:
%   spikePrs = setSpikePrs(mosaic)
%
% Description:
%    Converts the mosaic object into the appropriate format for the Pillow
%    spike-generating code.
% 
%    The properties of 'glmprs' from "testscript_GLM_coupled.m" found in
%    \pillow_code_GLM_v1_Feb2010, available at 
%    http://pillowlab.princeton.edu/code_GLM.html.
% 
%    Target output structure:
%      ggsim = 
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
%    mosaic   - Object. A rgc mosaic object.
%
% Outputs:
%    spikePrs - Struct. A structure of the format described above.
%
% Optional key/value pairs:
%    None.
%
% References:
%    * http://pillowlab.princeton.edu/code_GLM.html.

% History:
%    03/XX/16  JRG  (c) isetbio
%    05/30/19  JNM  Documentation pass

%% Get number of cells
nCells = size(mosaic.cellLocation);

%% Set linear temporal filter
% We do the linear temporal filtering before getting here, at the linear
% response level. We do not want any further temporal filtering, so we set
% the temporal responses to impulses.
spikePrs.k = zeros(length(mosaic.tCenter{1}), 1, 2);
spikePrs.k = zeros(20, nCells(1) * nCells(2), nCells(1) * nCells(2));
cellCtr = 0;
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
        spikePrs.k(20, cellCtr, cellCtr) = 1;
    end
end

%% Set DC value
spikePrs.dc = zeros(size(mosaic.sRFcenter));

%% Set post-spike filters
postSpikeFilter = mosaicGet(mosaic, 'postSpikeFilter');

% Reformat the postSpikeFilter into the key variable from the Pillow code
hlen = length(postSpikeFilter);
nCells = size(mosaic.cellLocation);
nCellsTotal = nCells(1)*nCells(2);
% spikeTimes = cell(nCells);
cellCtr = 0;
ih = zeros(nCellsTotal, nCellsTotal, hlen);
for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        cellCtr = cellCtr+1;
        ih(cellCtr, cellCtr, :) = postSpikeFilter;
    end
end
    
ih = permute(ih, [3 2 1]); % flip 2nd & 3rd dimensions

spikePrs.ih = ih;

%% Set time samples
spikePrs.iht = mosaic.dt*(1:length(spikePrs.ih));

%% Set interpolation
spikePrs.dt = mosaic.dt;

%% Set nonlinearity
spikePrs.nlfun = mosaic.generatorFunction;

end
