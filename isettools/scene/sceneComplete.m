function scene = sceneComplete(scene)
% Complete the scene parameters (geometry and spatial) after initialization
%
% Syntax:
%	scene = sceneComplete(scene);
%
% Description
%    This is the code at the end of sceneCreate, after the scene is
%    basically constructed. It initialize the Geomtry and Spatial
%    parameters. At some point, we might use this function rather than the
%    code in sceneCreate.
%
% Inputs:
%    scene - A scene structure
%
% Outputs:
%    scene - The modified scene structure
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    sceneInit, sceneCreate, sceneHarmonic
%

% History:
%    xx/xx/17  BW   ISETBIO Team, 2017
%    01/26/18  jnm  Formatting

%% Initialize scene geometry, spatial sampling if not already set above
scene = sceneInitGeometry(scene);
scene = sceneInitSpatial(scene);

% Scenes are initialized to a mean luminance of 100 cd/m2. The illuminant
% is adjusted so that dividing the radiance (in photons) by the illuminant
% (in photons) produces the appropriate peak reflectance (default = 1).
%
% Also, a best guess is made about one known reflectance.
% if isfield(scene, 'data')
if checkfields(scene, 'data', 'photons')
    if isempty(sceneGet(scene, 'knownReflectance')) && ...
            checkfields(scene, 'data', 'photons')       
        % If there is no illuminant yet, set the illuminant to equal
        % photons at 100 cd/m2
        wave = sceneGet(scene, 'wave');
        if isempty(sceneGet(scene, 'illuminant'))
            il = illuminantCreate('equal photons', wave, 100);
            scene = sceneSet(scene, 'illuminant', il);
        end
        
        % There is no knownReflectance, so we set the peak radiance to a
        % reflectance of 0.9.
        v = sceneGet(scene, 'peakRadianceAndWave');
        idxWave = find(wave == v(2));
        p = sceneGet(scene, 'photons', v(2));
        [~, ij] = max2(p);
        v = [0.9 ij(1) ij(2) idxWave];
        scene = sceneSet(scene, 'knownReflectance', v);
    end
    
    luminance = sceneCalculateLuminance(scene);
    scene = sceneSet(scene, 'luminance', luminance);
    
    % This routine also adjusts the illumination level to be consistent
    % with the reflectance and scene photons.
    scene = sceneAdjustLuminance(scene, 100);
    
    if ieSessionGet('gpu compute')
        p = sceneGet(scene, 'photons');
        scene = sceneSet(scene, 'photons', gpuArray(p));
    end
end
