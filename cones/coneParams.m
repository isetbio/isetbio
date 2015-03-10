function p = coneParams(coneAperture,scene)
%Return default cone parameters
%
%   p = coneParams(coneAperture,scene)
%
% Example:
%   scene = sceneCreate;
%   p = coneParams(2e-6,scene);
%   scene = sceneSet(scene,'fov',1);
%   p = coneParams(2e-6,scene);
%
%   p = coneParms;
%   h = sensorCreate('human',[],p);
%
% See also: sensorCreateConeMosaic, sensorCreate;
%
% (c) Stanford VISTA, Wandell, 2010

if notDefined('coneAperture'), coneAperture = 1.5e-6; end

% Cone aperture in meters.
% Central human cones are 1.5 um in the fovea, 3um in the periphery
p.coneAperture = [coneAperture coneAperture];

% Blank, L, M, S ratios
p.rgbDensities = [0 0.6, 0.3, 0.1];

% Could be mouse, I suppose
p.species = 'human'; 

% Random number
if exist('rng','file'), p.rSeed = rng;
else                    p.rSeed = randi([1,5000]);
end

if notDefined('scene')
    % Default size is small.  Spans about 3/4 a degree.
    p.sz = [72 88];
else
    % Match the array size to the scene, but a little smaller to allow for
    % eye movements.
    vFov = sceneGet(scene,'vfov');
    hFov = sceneGet(scene,'hfov');

    % Suppose there are 280 um/deg then 
    degPerCone = p.coneAperture/280e-6;
    p.sz = floor([(vFov/degPerCone(1)) (hFov/degPerCone(1))]);
    p.sz = round(p.sz * .9);

end

return;
