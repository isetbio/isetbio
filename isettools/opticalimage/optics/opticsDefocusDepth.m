function oDist = opticsDefocusDepth(defocus, optics, imgPlaneDist)
% Compute depth in object space to achieve a particular defocus (diopters).
%
% Syntax:
%   oDist = opticsDefocusDepth(defocus, optics, [imgPlaneDist])
%
% Description:
%    Given an optics structure and an image plot, we compute where the
%    object must be (in meters) from the lens to be set at one of the
%    defocus values sent in.
%
% Inputs:
%    defocus      - Vector. List of defocus levels in diopters
%    optics       - Struct. Optics structure
%    imgPlaneDist - (Optional) Numeric. The distance of the image plane
%                   behind the lens. Must be further than the focal length.
%                   Default is to retrieve the focal length.
%
% Outputs:
%    oDist        - Vector. Object distances to achieve the various
%                   specified defocus levels.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   opticsDepthDefocus, s3d_DepthSpacing
%

% History:
%    xx/xx/11       Copyright ImagEval Consultants, LLC, 2011.
%    03/09/18  jnm  Formatting
%    06/27/19  JNM  Formatting update

% Examples:
%{
    optics = opticsCreate;
    defocus = (-2:.1:-0.1);
    imgPlaneDist = opticsGet(optics, 'focal length', 'm');
    oDist = opticsDefocusDepth(defocus, optics, imgPlaneDist)
    figure;
    plot(defocus, oDist, '-o');
    xlabel('Defocus')
    ylabel('Obj Dist')

    % This is the inverse function that takes depth into defocus
   opticsDepthDefocus(oDist, optics, imgPlaneDist)
%}

if notDefined('defocus'), error('Defocus in diopters required'); end
if notDefined('optics'), error('Optics required'); end
if notDefined('imgPlaneDIst')
    imgPlaneDist = opticsGet(optics, 'focal length', 'm');
end

% Focal length needed for depth
fLength = opticsGet(optics, 'focal length', 'm');
if imgPlaneDist < fLength
    error('imgPlaneDist should be larger than the focal length')
end

% Invert the first one to get imgDist from defocus
imgDist = 1 ./ (defocus + (1 / imgPlaneDist));

% Invert
%   imgDist = (fLength * objDist) ./ (objDist - fLength);
% To get objDist from imgDist
oDist = imgDist * fLength ./ (imgDist - fLength);

end
