function glmprs = setGLMprs(mosaic, varargin)
% Convert mosaic object to a pillow spike-generating code format.
%
% Syntax:
%   glmprs = setGLMprs(mosaic, [varargin])
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
%        glmprs =
%               k: [20x1x2 double]
%              dc: [3 3]
%              ih: [598x2x2 double]
%             iht: [598x1 double]
%              dt: 0.0100
%           nlfun: @exp
%        ihbasprs: [1x1 struct]
%
% Inputs:
%    mosaic   - Object. A rgc mosaic object.
%
% Outputs:
%    glmprs   - Struct. A structure containing the GLM PRS information in
%               the format described above.
%
% Optional key/value pairs:
%    coupling - Boolean. A boolean for whether or not there is coupling.
%

% History:
%    03/XX/16  JRG  (c) isetbio
%    06/07/19  JNM  Documentation pass

%% Parse
p = inputParser;
p.addRequired('mosaic');
p.addParameter('coupling', true, @islogical);
p.parse(mosaic, varargin{:});
coupling = p.Results.coupling;

% Parameters shared whether coupled or not
glmprs.nlfun = mosaic.generatorFunction;

global RefreshRate

if coupling
    %% Get number of cells
    nCells = mosaic.get('mosaic samples');
    nCellsTotal = nCells(1) * nCells(2);

    %% Set linear temporal filter
    % Because of the construction of simGLM and simGLMcpl, we compute the
    % linear spatial and temporal responses before passing to simGLM. In
    % order to avoid any temporal filtering internal to simGLM, we set the
    % temporal responses to impulses. Taking the convolution of each
    % neuron's linear response with an impulse then results in the original
    % linear response.

    % The 20 is the default size for the impulse response
    glmprs.k = zeros(20, nCellsTotal, nCellsTotal);

    for ii = 1:nCellsTotal, glmprs.k(20, ii, ii) = 1; end

    %% Set DC value
    glmprs.dc = zeros(nCells);

    %% Set coupling filters
%     ihcpl = mosaicGet(mosaic, 'couplingFilter');
%
%     hlen = length(ihcpl{1, 1});
%     cellCtr = 0;
%     ih = zeros(nCellsTotal, nCellsTotal, hlen);
%     for xcell = 1:nCells(1)
%         for ycell = 1:nCells(2)
%             cellCtr = cellCtr+1;
%
%             if isa(mosaic, 'rgcPhys')
%                 % @JRG - Move out to EJ repository
%                 coupledCells = mosaic.couplingMatrix{xcell, ycell};
%                 ih(cellCtr, coupledCells, :) = ...
%                      (horzcat(mosaic.couplingFilter{1, cellCtr}{:}))';
%                 ih(cellCtr, cellCtr, :) = ...
%                      ((mosaic.postSpikeFilter{1, cellCtr}))';
%             else
%                 ih(cellCtr, :, :) = reshape(...
%                      mosaic.couplingFilter{cellCtr}, nCellsTotal, hlen);
%             end
%         end
%     end
%
%     ih = permute(ih, [3 2 1]); % flip 2nd & 3rd dimensions
%     glmprs.ih = ih;

    ihcpl = mosaicGet(mosaic, 'couplingFilter');

    hlen = length(ihcpl{1, 1});
    cellCtr = 0;
    ih = zeros(6, nCellsTotal, hlen);
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            cellCtr = cellCtr + 1;

            if isa(mosaic, 'rgcPhys')
                % @JRG - Move out to EJ repository
                coupledCells = mosaic.couplingMatrix{xcell, ycell};
                ih(cellCtr, coupledCells, :) = ...
                    (horzcat(mosaic.couplingFilter{1, cellCtr}{:}))';
                ih(cellCtr, cellCtr, :) = ...
                    ((mosaic.postSpikeFilter{1, cellCtr}))';
            else
%                 ih(cellCtr, :, :) = reshape(...
%                    mosaic.couplingFilter{cellCtr}, nCellsTotal, hlen);
                clear ii jj indNZ
                indNZ = find(mosaic.couplingFilter{cellCtr}(:, :, 1) ~= 0);
                [ii, jj] = ind2sub([mosaic.get('mosaic samples')], indNZ);

                for iiind = 1:length(ii)
                    ih(iiind, cellCtr, :) = mosaic.couplingFilter{...
                        cellCtr}(ii(iiind), jj(iiind), :);
                end
                ihind(cellCtr, 1:length(indNZ)) = indNZ;
            end
        end
    end

    ih = permute(ih, [3 2 1]); % flip 2nd & 3rd dimensions
    glmprs.ih = ih;
    glmprs.ihind = ihind;
    %% Set time samples
    if~(isfield(mosaic, 'dt'))
        glmprs.iht = (1 / RefreshRate) * (1:size(ih, 1));
        glmprs.dt = (1 / RefreshRate);
    else
        % Set interpolation
        glmprs.iht = mosaic.dt  *(1:size(ih, 1));
        glmprs.dt = mosaic.dt;
    end

%% Set nonlinearity
else
    % No coupling. Simpler.
    nCells = [1, 1];
    nCellsTotal = 1;

    %% Set linear temporal filter
    % Because of the construction of simGLM and simGLMcpl, we compute the
    % linear spatial and temporal responses before passing to simGLM. In
    % order to avoid any temporal filtering internal to simGLM, we set the
    % temporal responses to impulses. The convolution of each neuron's
    % linear response with an impulse then results in the (desired)
    % original linear response.
    glmprs.k = zeros(20, nCellsTotal, nCellsTotal);
    for ii = 1:nCellsTotal, glmprs.k(20, ii, ii) = 1; end

    %% Set DC value
    glmprs.dc = zeros(nCells);

    %% Set post spike filters
    ih = ((mosaic.postSpikeFilter));
    glmprs.ih = ih;    

    %% Set time samples
    if~(isfield(mosaic, 'dt'))
        glmprs.iht = (1 / RefreshRate) * (1:size(ih, 1));
        glmprs.dt = (1 / RefreshRate);
    else
        % Set interpolation
        glmprs.iht = mosaic.dt * (1:size(ih, 1));
        glmprs.dt = mosaic.dt;
    end
end

end
