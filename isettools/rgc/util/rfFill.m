function fillIndices = rfFill(obj)
% Fill in empty indices in the provided rgc mosaic
%
% Syntax:
%   fillIndices = rfFill(obj)
%
% Description:
%    Fill in empty indices in a retinal ganglion cell mosaic.
%
% Inputs:
%    obj         - Object. A retinal ganglion cell mosaic object.
%
% Outputs:
%    fillIndices - Cell. A cell array object indicating the filled indices.
%
% Optional key/value pairs:
%    None.
%

nCells = size(obj.cellLocation);

for xcell = 1:nCells(1)
    for ycell = 1:nCells(2)
        fillIndices{xcell, ycell}  = ...
            find(obj.sRFcenter{xcell, ycell} > ...
            obj.rfDiaMagnitude{xcell, ycell});
    end
end
