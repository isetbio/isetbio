function dMap = oiPadDepthMap(scene, invert)
% Pad the scene depth map into the oi depth map
%
% Syntax:
%   dMap = oiPadDepthMap(scene)
%
% Description:
%    It is also possible to extract the center part of the depth map,
%    without the padding.
%
% Inputs:
%    scene  - Struct. A scene structure.
%    invert - (Optional) Boolean. A boolean indicating whether or not to
%             invert the padding. Default is false.
%
% Outputs:
%    dMap   - Cell array. The cell array containing the depth map.
%
% Optional key/value pairs:
%    None.
%

% Check variables
if notDefined('scene'), error('Scene is required.'); end
if notDefined('invert'), invert = 0; end

if ~invert
    dMap = sceneGet(scene, 'depth map');
    sSize = sceneGet(scene, 'size');
    padSize = round(sSize / 8);
    padSize(3) = 0;
    padval = 0;
    direction = 'both';
    dMap = padarray(dMap, padSize, padval, direction);
    % figure;
    % imagesc(dMap)
else
    dMap = oiGet(oi, 'depth map');
    [r, c] = size(dMap);
    center = round(r / 2);
    cDepth = dMap(center, :);
    hData = (cDepth > 0);
    hData = double(hData);
    center = round(c / 2);
    cDepth = dMap(:, center);
    vData = (cDepth > 0);
    vData = double(vData);
    gData = logical(hData(:) * vData(:)');
    figure;
    imagesc(gData)
    % oMap = dMap(gData);
    % figure;
    % imagesc(oMap)
end

end