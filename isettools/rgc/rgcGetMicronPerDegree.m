function micronPerDegree = rgcGetMicronPerDegree(scene,oi)
% Calculate the number of microns per degree of visual field
%
%  micronPerDegree = rgcGetMicronPerDegree(scene,oi)
%
% (c) Stanford Synapse Team 2010

f = oiGet(oi, 'optics flength');
d = sceneGet(scene,'distance');

% Thin lens formula
meterPerDegree = (1/(1/f - 1/d))*tand(1);

% converting to microns
micronPerDegree = meterPerDegree*1e6;

end