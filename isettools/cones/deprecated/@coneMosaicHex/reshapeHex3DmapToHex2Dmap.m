function hex2Dmap = reshapeHex3DmapToHex2Dmap(obj, hex3Dmap)
% Reshape the existing 3D Hex Map to a 2D Hex Map
%
% Syntax:
%   reshapeHex3DmapToHex2Dmap(obj, hex3Dmap)
%
% Description:
%    Take a 3D Hex Map and map it to a 2D Hex Map
%
% Inputs:
%    obj      - The cone mosaic hex object
%    hex3Dmap - The 3D Hex Map
%
% Outputs:
%    hex2Dmap - The created 2D Hex Map
%
% Optional key/value pairs:
%    None.
%
nonNullCones = obj.pattern(obj.pattern > 1);

iLdest = nonNullCones == 2;
iMdest = nonNullCones == 3;
iSdest = nonNullCones == 4;

iLsource = obj.pattern == 2;
iMsource = obj.pattern == 3;
iSsource = obj.pattern == 4;

timePointsNum = size(hex3Dmap, 3);
hex3Dmap = reshape(hex3Dmap, ...
    [size(obj.pattern, 1) * size(obj.pattern, 2) timePointsNum]);
hex2Dmap = zeros(numel(nonNullCones), timePointsNum, class(hex3Dmap));
hex2Dmap(iLdest, :) = hex3Dmap(iLsource, :);
hex2Dmap(iMdest, :) = hex3Dmap(iMsource, :);
hex2Dmap(iSdest, :) = hex3Dmap(iSsource, :);

end
