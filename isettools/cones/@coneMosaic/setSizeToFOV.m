function obj = setSizeToFOV(obj, fov, varargin)
% set size to fov
% Updates the cone mosaic size to match the FOV
%
%    cm.setSizeToFOV(fov,varargin)
%
% Inputs:
%   fov - 2-element vector for desired horizontal/vertical
%         field of view in degrees
%
% Parameters
%     sceneDist:    (Inf,   meters)
%     focalLength:  (0.017, meters)
%
% HJ ISETBIO Team 2016

% parse input
p = inputParser;
p.addRequired('fov', @isnumeric);
p.addParameter('sceneDist', inf, @isnumeric);
p.addParameter('focalLength', 0.017, @isnumeric);

p.parse(fov, varargin{:});
sceneDist = p.Results.sceneDist;
focalLength = p.Results.focalLength;
if isscalar(fov), fov = [fov fov]; end
fov = fov([2 1]); % flip the order to be vertical/horizontal

% compute current field of view
imageDist = 1/(1/focalLength - 1/sceneDist);
curFov = 2*atand([obj.height obj.width]/imageDist/2);

% set new size to object
obj.mosaicSize = ceil(obj.mosaicSize .* fov./curFov);
end