function glmprs = setGLMprs(mosaic,varargin)
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

%% Parse
p = inputParser;
p.addRequired('mosaic');
p.addParameter('coupling',true,@islogical);
p.parse(mosaic,varargin{:});
coupling = p.Results.coupling;

% Parameters shared whether coupled or not
glmprs.nlfun = mosaic.generatorFunction;

if coupling
    %% Get number of cells
    nCells = mosaic.get('mosaic size');
    nCellsTotal = nCells(1)*nCells(2);
    
    %% Set linear temporal filter
    % Because of the construction of simGLM and simGLMcpl, we compute the
    % linear spatial and temporal responses before passing to simGLM. In
    % order to not do any temporal filtering internal to simGLM, we set the
    % temporal responses to impulses. The convolution of each neuron's
    % linear response with an impulse results in the original linear
    % response.
    
    % The 20 is the default size for the impulse response
    glmprs.k = zeros(20,nCellsTotal,nCellsTotal);
    
    for ii=1:nCellsTotal
        glmprs.k(20,ii,ii) = 1;
    end
    
    
    %% Set DC value
    glmprs.dc = zeros(nCells);
    
    %% Set coupling filters
    ihcpl = mosaicGet(mosaic, 'couplingFilter');
    
    hlen = length(ihcpl{1,1});
    cellCtr = 0;
    ih = zeros(nCellsTotal, nCellsTotal, hlen);
    for xcell = 1:nCells(1)
        for ycell = 1:nCells(2)
            cellCtr = cellCtr+1;
            
            if isa(mosaic, 'rgcPhys')
                % @JRG - Move out to EJ repository
                coupledCells = mosaic.couplingMatrix{xcell,ycell};
                ih(cellCtr,coupledCells,:) = (horzcat(mosaic.couplingFilter{1,cellCtr}{:}))';
                ih(cellCtr,cellCtr,:) = ((mosaic.postSpikeFilter{1,cellCtr}))';
            else
                ih(cellCtr,:,:) = reshape(mosaic.couplingFilter{cellCtr}, nCellsTotal,hlen);
            end
        end
    end
    
    ih = permute(ih,[3 2 1]); % flip 2nd & 3rd dimensions
    glmprs.ih = ih;
    
    %% Set time samples
    if~(isfield(mosaic,'dt'))
        glmprs.iht = .01*(1:size(ih,1));
        glmprs.dt = .01;
    else
        % Set interpolation
        glmprs.iht = mosaic.dt*(1:size(ih,1));
        glmprs.dt = mosaic.dt;
    end
    
    %% Set nonlinearity
else
    % No coupling.  Simpler.
    nCells = [1,1];
    nCellsTotal = 1;
    
    %% Set linear temporal filter
    % Because of the construction of simGLM and simGLMcpl, we compute the
    % linear spatial and temporal responses before passing to simGLM. In
    % order to not do any temporal filtering internal to simGLM, we set the
    % temporal responses to impulses. The convolution of each neuron's
    % linear response with an impulse results in the original linear
    % response. 
    
    glmprs.k = zeros(20,nCellsTotal,nCellsTotal);
    
    for ii=1:nCellsTotal
        glmprs.k(20,ii,ii) = 1;
    end
    
    
    %% Set DC value
    glmprs.dc = zeros(nCells);
    
    %% Set post spike filters

    ih = ((mosaic.postSpikeFilter));
    glmprs.ih = ih;    
    
    %% Set time samples
    if~(isfield(mosaic,'dt'))
        glmprs.iht = .01*(1:size(ih,1));
        glmprs.dt = .01;
    else
        % Set interpolation
        glmprs.iht = mosaic.dt*(1:size(ih,1));
        glmprs.dt = mosaic.dt;
    end
end

end

