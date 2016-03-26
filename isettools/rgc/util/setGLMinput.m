function glminput = setGLMinput(response)
% Reformats a mosaic's linear response into the simGLM or simGLMcpl format
%
% This is used for running Pillow's code.  At the moment we use it for
% other models, but that should be eliminated over time.
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
glminput = zeros(nCellsTotal,length(response{1,1,1}));

for i = 1:nCells(1)
    for j = 1:nCells(2)        
        cellCtr = cellCtr+1;
        %  glminput(cellCtr,:) = (squeeze(sum(sum(response{i,j}(:,:,:,1),1),2)))';
        glminput(cellCtr,:) = response{i,j,1};
    end
end

end