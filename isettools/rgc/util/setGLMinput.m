function glminput = setGLMinput(response)

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
