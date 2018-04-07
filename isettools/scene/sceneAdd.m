function scene1 = sceneAdd(scene1, scene2, addFlag)
% Add together the photons from two scenes
%
% Syntax:
%   scene = sceneAdd(scene1, scene2, [addFlag])
%
% Description:
%    Add two the radiance from two scenes that match in all ways (row, col, 
%    wavelength, so forth).
%
% Inputs:
%    scene1  - Original scene
%    scene2  - Scene to add to the original scene.
%    addFlag - (Optional) String indicating whether to add the two scenes,
%              or remove the spatial mean first. Default 'add'. Options for
%              the argument are:
%       add               - Add the radiance data from the scenes. Default.
%       removespatialmean - Remove mean from scene2 and then add scenes.
%
% Outputs:
%    scene1  - The combined scene.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/12       Copyright Imageval 2012
%    12/28/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match Wiki.


%% We should do more parameter checking, I suppose
if notDefined('scene1'), error('scene 1 required'); end
if notDefined('scene2'), error('scene 2 required'); end
if notDefined('addFlag'), addFlag = 'add'; end


%% Get the photons and do the right thing
p = sceneGet(scene1, 'photons');
s = sceneGet(scene2, 'photons');
nWave = sceneGet(scene2, 'nwave');

addFlag = ieParamFormat(addFlag);
switch addFlag
    case 'add'
    case 'removespatialmean'
        for ii = 1:nWave
            s(:, :, ii) = s(:, :, ii) - mean(mean(s(:, :, ii)));
        end
    otherwise
end

scene1 = sceneSet(scene1, 'photons', p + s);

end