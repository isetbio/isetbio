function [scene, peakReflectance] = sceneAdjustReflectance(scene, prct)
% Scale the data so that the peak reflectance is 1
%
% Syntax:
%   [scene, peakReflectance] = sceneAdjustReflectance(scene, [prct])
%
% Description:
%    We scale the illuminant to a higher level to adjust peak reflectance.
%
%    By default we find the highest reflectance (percentile = 100). If this
%    value is < 1, we do nothing. If the value exceeds one, we scale the
%    illuminant so that the max reflectance is 1.
%
%    If there are specularities and this scaling makes the reflectances too
%    small. To allow for these cases, choose a percentile less than 100, 
%    say 98.
%
%    N.B. The source contains executable examples of usage, which can be
%    accessed by typing 'edit sceneAdjustReflectance.m' into MATLAB's
%    command window.
%
%    There are examples in the code. Type 'edit sceneAdjustReflectance'
%    into the Command Window to access.
%
% Inputs:
%    scene           - A scene structure.
%    prct            - (Optional) The percentile which we scale the
%                      reflectance by. Default is 100.
%
% Outputs:
%    scene           - The modified scene structure.
%    peakReflectance - The peak reflectance of the scene structure.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    sceneAdjustIlluminant, sceneAdjustLuminance
%

% History:
%    xx/xx/12       (c) Imageval, 2012
%    12/29/18  jnm  Formatting
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    fullFileName = fullfile(isetbioDataPath, 'images', ...
        'multispectral', 'StuffedAnimals_tungsten-hdrs.mat');
    scene = sceneFromFile(fullFileName, 'multispectral');
    scene = sceneAdjustReflectance(scene);
    vcAddAndSelectObject('scene', scene);
    sceneWindow;

    scene = sceneAdjustReflectance(scene, 98);
    vcAddAndSelectObject('scene', scene);
    sceneWindow;
%}

if notDefined('scene'), error('Scene required.'); end
if notDefined('prct'), prct = 100; end

% Find the current peak reflectance value
r = sceneGet(scene, 'reflectance');

% Find the value that we want to be at unit reflectance. Because of
% specularities, it may be that some locations will have > 1 reflectance.
if prct == 100
    peakReflectance = max(r(:));
else
    % Faster band-wise. Not quite the same as the global percentile. This
    % is the band-wise percentile.
    nWave = sceneGet(scene, 'nwave');
    peakReflectanceBand = zeros(nWave, 1);
    str = sprintf('Finding %d percentile reflectance', prct);
    showBar = ieSessionGet('waitbar');
    if showBar, h = waitbar(0, str); end
    for ii = 1:nWave
        if showBar, waitbar(ii / nWave, h); end
        tmp = r(:, :, ii);
        peakReflectanceBand(ii) = iePrctile(tmp(:), prct);
    end
    if showBar, close(h); end
    peakReflectance = iePrctile(peakReflectanceBand, prct);
end

if peakReflectance > 1
	newIlluminant = sceneGet(scene, 'illuminantPhotons') * peakReflectance;
    scene = sceneSet(scene, 'illuminantPhotons', newIlluminant);
else
    disp('All reflectance values are < 1');
end

end