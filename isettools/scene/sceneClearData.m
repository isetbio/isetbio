function scene = sceneClearData(scene)
%  Clear the scene data entries. 
%
%   scene = sceneClearData(scene)
%
% Copyright ImagEval Consultants, LLC, 2003.

bitDepth = sceneGet(scene,'bit depth');

scene = sceneSet(scene,'data',[]);
scene = sceneSet(scene,'bit depth',bitDepth);

%% Illuminant
% scene = sceneSet(scene,'illuminant energy',[]);

end