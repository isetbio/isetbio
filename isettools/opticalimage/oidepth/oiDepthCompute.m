function [oiD, D] = oiDepthCompute(oi, scene, imageDist, depthEdges, ...
    cAberration, displayFlag)
% Compute cell array of OI for different depth planes from scene
%
% Syntax:
%   [oiD, D] = oiDepthCompute(oi, scene, imageDist, depthEdges, ...
%       cAberration, displayFlag)
%
% Description:
%    Compute the cell array of OI for different depth planes from a scene.
%
% Inputs:
%    oi          - Struct. An optical image structure.
%    scene       - Struct. A scene structure.
%    imageDist   - (Optional) Numeric. The distance of the image plane
%                  behind the lens. Default is the focal length.
%    depthEdges  - Vector. A vector containing distances from the lens. The
%                  number of depthEdges defines the number of OIs that will
%                  be computed.
%    cAberration - (Optional) Numeric. Defocus in diopters (chromatic
%                  aberration, longitudinal) for each wavelength. Default
%                  is [].
%    displayFlag - Boolean. Display the defocused images is shown in the
%                  oiWindow. Default is true.
%
% Outputs:
%    oiD         - Cell. A cell array of irradiance images. There is one
%                  image per distance (m) in depthEdges. The defocused
%                  image is calculated for each distance. This cell array
%                  of irradiance images is combined (oiDepthCombine) into a
%                  single image. The combination is based on picking out
%                  the pixels at the appropriate depth.
%    D           - Vector. The defocus per each image in oiD.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - s3dRenderDepthDefocus return specified in that function
%      is oi, oiD, D ... however, we are taking the first and third of
%      those arguments and calling them oiD and D. Can we double check that
%      we are parsing this properly?]
%
% See Also:
%   oiDepthCombine, oiDepthSegmentMap, s_opticsDepthScene
%

% History:
%    xx/xx/11       Copyright ImagEval Consultants, LLC, 2011.
%    03/28/18  jnm  Formatting
%    07/01/19  JNM  Formatting update

if notDefined('oi'), error('oi required'); end
if notDefined('scene'), error('scene required'); end
if notDefined('imageDist')
    imageDist = opticsGet(oiGet(oi, 'optics'), 'focal length', 'm');
end
if notDefined('depthEdges'), error('depthEdges required'); end
if notDefined('cAberration'), cAberration = []; end
if notDefined('displayFlag'), displayFlag = 1; end

% Set the scene map to a single depth. We  sweep through the depthEdges for
% the whole scene.
oMap  = sceneGet(scene, 'depth map');
[r, c] = size(oMap);
dMap  = ones(r, c);

% Cell array of oi images at different depths
oiD = cell(1, length(depthEdges));

% Loop and create the OI structures
for ii = 1:length(depthEdges)
    scene = sceneSet(scene, 'depth map', dMap * depthEdges(ii));
    [oiD{ii}, ~, D] = s3dRenderDepthDefocus(...
        scene, oi, imageDist, [], cAberration);
    oiD{ii} = oiSet(oiD{ii}, 'name', sprintf('Defocus %.2f', D));
    if displayFlag
        vcAddAndSelectObject(oiD{ii});
        oiWindow
    end
end

end