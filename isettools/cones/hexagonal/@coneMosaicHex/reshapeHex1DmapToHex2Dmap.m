function hex2Dmap = reshapeHex1DmapToHex2Dmap(obj, hex1Dmap)
% Reshape the existing 1D Hex Map to a 2D Hex Map
%
% Syntax:
%   reshapeHex1DmapToHex2Dmap(obj, hex1Dmap)
%
% Description:
%    Take a 1D Hex Map returned by computeForOISequence and map it to a 
%    2D Hex Map
%
% Inputs:
%    obj      - The cone mosaic hex object
%    hex1Dmap - The 1D Hex Map
%
% Outputs:
%    hex2Dmap - The created 2D Hex Map
%
% Optional key/value pairs:
%    None.
%

hex2Dmap = 0 * obj.pattern;

nonNullCones = obj.pattern(obj.pattern > 1);
iLsource = nonNullCones == 2;
iMsource = nonNullCones == 3;
iSsource = nonNullCones == 4;

iLdest = obj.pattern == 2;
iMdest = obj.pattern == 3;
iSdest = obj.pattern == 4;

hex2Dmap(iLdest) = hex1Dmap(iLsource);
hex2Dmap(iMdest) = hex1Dmap(iMsource);
hex2Dmap(iSdest) = hex1Dmap(iSsource);

end

