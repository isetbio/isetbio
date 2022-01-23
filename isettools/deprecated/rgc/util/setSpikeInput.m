function spikeInput = setSpikeInput(response)
% Reformats a mosaic's linear response into the simGLM or simGLMcpl format
%
% Syntax:
%   spikeInput = setSpikeInput(response)
%
% Description:
%    Reformat a mosaic's linear response into either the simGLM or the
%    simGLMcpl format.
%
%    This is used for running Pillow's code.  At the moment we use it for
%    other models, but that should be eliminated over time.
%
% Inputs:
%    response   - <Type>. A rgc mosaic linear response.
%
% Outputs:
%    spikeInput - <Type>. A simGLM or simGLMcpl formatted response.
%
% Optional key/value pairs:
%    None.
%
% References:
%    * \pillow_code_GLM_v1_Feb2010, available at 
%      http://pillowlab.princeton.edu/code_GLM.html.
%

% History:
%    03/XX/16  JRG  (c) isetbio
%    05/31/19  JNM  Documentation pass

%% Set glminput to linear response
nCells = size(response);
nCellsTotal = nCells(1) * nCells(2);
cellCtr = 0;
spikeInput = zeros(nCellsTotal, length(response{1, 1, 1}));

for i = 1:nCells(1)
    for j = 1:nCells(2)        
        cellCtr = cellCtr + 1;
        % glminput(cellCtr, :) = ...
        %     (squeeze(sum(sum(response{i, j}(:, :, :, 1), 1), 2)))';
        spikeInput(cellCtr, :) = response{i, j, 1};
    end
end

end