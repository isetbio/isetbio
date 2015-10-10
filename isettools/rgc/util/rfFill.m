function fillIndices = rfFill(obj)
nCells = size(obj.cellLocation);

for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        fillIndices{xcell,ycell}  = find(obj.sRFcenter{xcell,ycell} > obj.rfDiaMagnitude{xcell,ycell});
    end
end


       

