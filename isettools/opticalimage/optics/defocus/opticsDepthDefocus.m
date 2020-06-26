function [D, imgDist] = opticsDepthDefocus(objDist, optics, imgPlaneDist)
% Compute defocus in diopters for vector of object distances (meters)
%
% Syntax:
%   [D, imgDist] = opticsDepthDefocus(objDist, [optics], [imgPlaneDist])
%
% Description:
%    The lensmaker's equation specifies the relationship between object and
%    image plane distances given the focal length of a thin lens. The
%    lensmaker's equation for a thin lens is
%
%       1 / objDist + 1 / imgDist = 1 / focalLength
%
%    See the thin lens equation description from Wikipedia.
%    http://en.wikipedia.org/wiki/Lens_(optics)
%
%    The equation can be used for various purposes, including specifying
%    the defocus of an object in various imaging conditions.
%
%    For example,
%      * Objects at infinity are imaged in the focal plane.
%      * If the image plane is at the focal length, then closer objects
%        will be in defocus and we can assess the degree of defocus. These
%        will all be positive.
%      * If the image plane distance differs from the focal length, we can
%        compute the defocus as a function of object distance. Object
%        distances closer and further than the in-focus distance will have
%        negative and positive defocus.
%
%    There are additional lensmaker's formulae for general lenses with
%    specified curvature and index of refraction measurements. There are
%    also depth of field formulae in Wikipedia and elsewhere.
%
%    Simple alegebraic observations about the lensmaker's equation
%
%    It is convenient to express the object distance as a multiple of focal
%    lengths, objDist = N * fLength. In that case
%
%      1/imgDist = 1 / fLength - 1 / (N * fLength)
%      1/imgDist = (1 / fLength) * (1 - 1 / N)
%      imgDist   = fLength / (1 - (1 / N))
%      imgDist   = fLength / ((N - 1) / N)
%      imgDist   = fLength * (N / (N - 1))
%
%    This expresses how the imgDist approaches the focal length as N gets
%    large. The difference between the image plane and the focal plane is
%
%      errDist = fLength - fLength * (N / (N - 1))
%      errDist = fLength * (1 - (N / (N - 1))
%
%    The blurring caused by this difference in best image plane depends on
%    the pupil aperture as well as the focal length. For small pupil
%    apertures, there is less of a penalty in defocus (e.g., a pinhole).
%    The significance of the pupil aperture is captured in other
%
%    This function contains examples of usage inline. To access these, type
%    'edit opticsDepthDefocus.m' into the Command Window.
%
% Inputs:
%    objDist      - Vector. Vector of object distances in meters.
%    optics       - (Optional) Struct. Optics structure. Default uses
%                   vcGetObject to pull an existing optics.
%    imgPlaneDist - (Optional) Numeric. Image plane distance (meters).
%                   Default is to use the focal length.
%
% Outputs:
%    D            - Vector. Defocus in diopters
%    imgDist      - Vector. The image distance at which the obj at objDist
%                   is in focus.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   opticsDefocusedMTF, defocusMTF, humanCore
%

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010.
%    03/16/18  jnm  Formatting
%    06/28/19  JNM  Documentation update

% Examples:
%{
    optics = opticsCreate;
    fLength = opticsGet(optics, 'focal length', 'm');
    nSteps = 500;
    objDist = linspace(fLength * 1.5, 50 * fLength, nSteps);
    D0 = opticsGet(optics, 'power');

    % Show defocus relative to total lens power
    figure(1)
    [D, imgDist] = opticsDepthDefocus(objDist, optics);
    plot(objDist, D / D0);
    xlabel('Distance to object (m)');
    ylabel('Relative dioptric error')

    % Object and image distances
    figure(2)
    [D, imgDist] = opticsDepthDefocus(objDist, optics);
    plot(objDist, imgDist);
    xlabel('Distance to object (m)');
    ylabel('Image dist (m)')
    line([objDist(1) objDist(end)], [fLength fLength], ...
        'Color', 'k', 'linestyle', '--')

    % Defocus with respect to image plane different from focal image plane
    figure(3)
    [D, imgDist] = opticsDepthDefocus(objDist, optics, 2 * fLength);
    plot(objDist / fLength, D);
    xlabel('Distance to object (re: Focal length)');
    ylabel('Dioptric error (1/m)')
    [v, ii] = min(abs(D));
    fprintf('In focus objDist:  %.3f re fLength %.3f\n', ...
        objDist(ii), objDist(ii) / fLength);
%}

if notDefined('objDist'), error('No distance specified'); end
if notDefined('optics'), optics = vcGetObject('optics'); end
fLength = opticsGet(optics, 'focal length', 'm');

if notDefined('imgPlaneDist'), imgPlaneDist = fLength; end
if imgPlaneDist < fLength
    error('Virtual image: img plane closer than fLength - not computed');
end

% Compute the image distance for various object distances

% Lensmaker's equation for a thin lens can be written:
%    1 / imgDist = (1 / fLength - 1 ./ objDist);
% so
%    imgDist = (1 / fLength - 1 ./ objDist) .^ -1;
% and ....
imgDist = (fLength * objDist) ./ (objDist - fLength);

% Compute the defocus - this is the dioptric power of a lens that would
% shift the image from its current distance (imgDist) to the desired image
% plane (imgPlaneDist).
D = (1 ./ imgDist) - (1 / imgPlaneDist);

%  figure(1);
%  plot(objDist, D)
%  xlabel('Distance to object (m)');
%  ylabel('Dioptric error (1/m)')
%  figure(1)

end