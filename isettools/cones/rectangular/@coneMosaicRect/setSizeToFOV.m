function obj = setSizeToFOV(obj, fov, varargin)
% Updates the cone mosaic size to match the FOV
%
% Syntax:
%   setSizeToFOV(obj, fov, [varargin])
%
% Description:
%    The functoin updates the cone mosaic size to match the field of view.
%    All input distances are assumed to be in meters.
%
% Inputs:
%    obj            - The cone mosaic object to modify
%    fov            - 2-element vector for desired horizontal/vertical
%                     field of view in degrees. If a scalar is passed, fov
%                     is taken to be square.
%
% Outputs:
%    obj            - The updated cone mosaic object
%
% Optional key/value pairs:
%     'sceneDist'   - Distance to scene in m. (default, Inf))
%     'focalLength' - Focal length assumed for eye in m. (default, 0.017)
%
% Notes:
%    * [Note: DHB - CAN THE DEFAULT FOCAL LENGTH FOR THE EYE BE OBTAINED
%      FROM THE CONEMOSAIC OBJECT?]
%    * [Note: DHB - IS INF A GOOD CHOICE FOR THE DEFAULT SCENE DISTANCE?
%      EXPLAIN IMPLICATIONS OF THIS CHOICE.]
%

% History:
%    xx/xx/16  HJ   ISETBIO Team 2016
%    02/19/18  jnm  Formatting

% parse input
p = inputParser;
p.addRequired('fov', @isnumeric);
p.addParameter('sceneDist', Inf, @isnumeric);
p.addParameter('focalLength', 0.017, @isnumeric);

p.parse(fov, varargin{:});
sceneDist = p.Results.sceneDist;
focalLength = p.Results.focalLength;

% Handle case where a scalar fov is passed. Interpret as a square fov.
if isscalar(fov), fov = [fov fov]; end

% This routine thinks in vertical/horizontal, flip input so code below
% works.
fov = fov([2 1]); 

% Compute current field of view
imageDist = 1 / (1 / focalLength - 1 / sceneDist);
curFov = 2 * atand([obj.height obj.width] / imageDist / 2);

% Set new size to object
obj.mosaicSize = ceil(obj.mosaicSize .* fov ./ curFov);

end