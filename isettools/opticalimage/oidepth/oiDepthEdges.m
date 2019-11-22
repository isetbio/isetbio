function [depthEdges, imageDist, oDefocus] = ...
    oiDepthEdges(oi, defocus, inFocusDepth)
% Determine depth edges to achieve defocus values
%
% Syntax:
%   [depthEdges, imageDist, oDefocus] = ...
%       oiDepthEdges(oi, defocus, inFocusDepth)
%
% Description:
%    Determine the depth edges to achieve the desired defocus values.
%
%    The defocus range should always start out as negative because when the
%    image plane is at the focal length infinite distance is in focus and
%    everything closer has negative defocus in diopters.
%
%    Some example values for an oiDepthEdges call:
%       defocus = -1.2:.2:0;   % Defocus when at image in focal plane
%       inFocusDepth = 5;      % Desired in focus depth (m)
%
% Inputs:
%    oi           - Struct. ISET optics. Optical image structure.
%    defocus      - Vector. Vector of defocus values.
%    inFocusDepth - Numeric. Depth for perfect focus.
%
% Outputs:
%    depthEdges   - Depths that achieve the relative defocus spacing when
%                   image is in the focal plane.
%    imageDist    - Image distance for the optics to achieve a best focus
%                   at inFocusDepth.
%    oDefocus     - The defocus at these depthEdges when the image plane is
%                   at imageDist.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/11       Copyright ImagEval Consultants, LLC, 2011.
%    03/28/18  jnm  Formatting

% Examples:
%{
    oi = oiCreate;
    optics = oiGet(oi, 'optics');
    f = opticsGet(optics, 'focal length', 'm');
    optics = opticsSet(optics, 'focal length', 5 * f);

    [depEdge, imDist, oDefocus] = oiDepthEdges(oi, [-1.2:.2:0], 5)
%}

optics = oiGet(oi, 'optics');
fLength = opticsGet(optics, 'focal length');
defocus(defocus >= 0) = -0.01;

% Calculate depths assuming image plane at focal length.
depthEdges = opticsDefocusDepth(defocus, optics, fLength);
[~, idx] = min(abs(inFocusDepth - depthEdges));

% If the user wants a particular image depth in focus, tell them where the
% image plane should be.
oDist = depthEdges(idx);
[~, imageDist] = opticsDepthDefocus(oDist, optics, fLength);

% This is the defocus for each depth when the image plane is at imageDist.
oDefocus = opticsDepthDefocus(depthEdges, optics, imageDist);

end