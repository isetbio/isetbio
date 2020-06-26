function idx = oiDepthSegmentMap(oiDmap, depthEdges)
% Depth segment map
%
% Syntax:
%   idx = oiDepthSegmentMap(oiDmap, depthEdges)
%
% Description:
%    Create the depth segment map for the provided oi depth map.
%
% Inputs:
%    oiDmap     - Cell Array. A cell array of OIs at different depths.
%    depthEdges - Vector. Depths that achieve the relative defocus spacing
%                 when image is in the focal plane.
%
% Outputs:
%    idx        - Vector. Indexes of map segments.
%
% Optional key/value pairs:
%    None.
%

nEdges = length(depthEdges);
[r, c] = size(oiDmap);
vMap = zeros(r, c, nEdges);

for ii = 1:nEdges, vMap(:, :, ii) = oiDmap - depthEdges(ii); end
% imagesc(oiDmap)

[v, idx] = min(abs(vMap), [], 3);
% figure;
% imagesc(idx)

return
