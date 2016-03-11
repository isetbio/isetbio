function glminput = setGLMinput(response)
% Converts the mosaic object's linear response into a format that can be
% passed to simGLM or simGLMcpl.
% 
% % \pillow_code_GLM_v1_Feb2010, available at 
% http://pillowlab.princeton.edu/code_GLM.html.
% 
% 
% 3/2016 JRG (c) isetbio

%% Set glminput to linear response
nCells = size(response);
nCellsTotal = nCells(1)*nCells(2);
cellCtr = 0;

for i = 1:nCells(1)
    for j = 1:nCells(2)        
        cellCtr = cellCtr+1;
%         glminput(cellCtr,:) = (squeeze(sum(sum(response{i,j}(:,:,:,1),1),2)))';
        glminput(cellCtr,:) = response{i,j,1};
    end
end
