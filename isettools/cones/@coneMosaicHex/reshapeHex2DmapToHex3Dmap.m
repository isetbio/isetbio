function hex3Dmap = reshapeHex2DmapToHex3Dmap(obj, hex2Dmap)
% Reshape the existing 2D Hex Map to a 3D Hex Map
%
% Syntax:
%   reshapeHex2DmapToHex3Dmap(obj, hex2Dmap)
%
% Description:
%    Take a 2D Hex Map and map it to a 3D Hex Map
%
% Inputs:
%    obj      - The cone mosaic hex object
%    hex2Dmap - The 2D Hex Map
%
% Outputs:
%    hex3Dmap - The created 3D Hex Map
%
% Optional key/value pairs:
%    None.
%
nonNullCones = obj.pattern(obj.pattern > 1);
iLsource = nonNullCones == 2;
iMsource = nonNullCones == 3;
iSsource = nonNullCones == 4;

iLdest = obj.pattern == 2;
iMdest = obj.pattern == 3;
iSdest = obj.pattern == 4;

timePointsNum = size(hex2Dmap, 2);
hex3Dmap = zeros(size(obj.pattern, 1) * size(obj.pattern, 2), ...
    timePointsNum, class(hex2Dmap));
hex3Dmap(iLdest, :) = hex2Dmap(iLsource, :);
hex3Dmap(iMdest, :) = hex2Dmap(iMsource, :);
hex3Dmap(iSdest, :) = hex2Dmap(iSsource, :);
hex3Dmap = reshape(hex3Dmap, size(obj.pattern, 1), ...
    size(obj.pattern, 2), timePointsNum);
end
