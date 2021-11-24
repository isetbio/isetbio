function coneLocsHexGrid = coneLocsHexGridAlignedWithSerializedConeMosaicResponse(obj)
% Ordered coneLocsHexGrid corresponds to the serialized cone response vector
%
% Syntax:
%    coneLocsHexGrid = coneLocsHexGridAlignedWithSerializedConeMosaicResponse(obj)
%
% Description:
%    Order the coneLocsHexGrid so it corresponds to the serialized cone
%    response vector
%
% Inputs:
%    obj - The cone mosaic hex object
%
% Outputs:
%    coneLocsHexGrid - the serialized coneLocsHexGrid
%
% Optional key/value pairs:
%   None

% History:
%    06/07/19  NPC  ISETBIO TEAM, 2015

    x = obj.patternSupport(:,:,1);
    y = obj.patternSupport(:,:,2);
    x = obj.reshapeHex3DmapToHex2Dmap(x);
    y = obj.reshapeHex3DmapToHex2Dmap(y);
    coneLocsHexGrid(:,1) = x;
    coneLocsHexGrid(:,2) = y;
end
