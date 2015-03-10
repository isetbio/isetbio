function fovS = getNewFov(fov,oiType,scene,oi)
% Computes the new fov necessary to use to compute a patch without being
% too much affect at the border.
% 
% fovS = getNewFov(fov,oiType,[scene],[oi]);
% 
% (c) Synapse Stanford Team - 2010

% Compute the new fov using the fact that the RF is limited at 20 um,
% reversing the blurring of the lens is done using an empirical function.

%% no argument cases:
if notDefined('scene')
    scene.distance = 1.2; % default scene distance in iset
end
if notDefined('oi')
    oi = rgcOiCreate(oiType);
end
RFSize = 20; % um

micronPerDegree = rgcGetMicronPerDegree(scene,oi);

fovp = fov + 2*(RFSize-1)/micronPerDegree;

% added 5 um for safety on each side
fovS = fovp+2*5/micronPerDegree;