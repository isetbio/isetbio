function scene = sceneClearData(scene)
%  Clear the scene data entries. 
%
% Syntax:
%   scene = sceneClearData(scene)
%
% Description:
%    Clear the scene data entries
%
% Inputs:
%    scene - The scene structure
%
% Outputs:
%    scene - The modified scene structure
%
% Optional key/value pairs:
%    None.

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    01/10/18  jnm  Formatting

bitDepth = sceneGet(scene, 'bit depth');

scene = sceneSet(scene, 'data', []);
scene = sceneSet(scene, 'bit depth', bitDepth);

%% Illuminant
% scene = sceneSet(scene, 'illuminant energy', []);

end