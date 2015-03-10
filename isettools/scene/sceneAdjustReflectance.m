function [scene, peakReflectance] = sceneAdjustReflectance(scene,prct)
%Scale the data so that the peak reflectance is 1
%
%   [scene, peakReflectance] = sceneAdjustReflectance(scene,[prct = 100])
%
% We scale the illuminant to a higher level to adjust the peak reflectance.
%
% By default we find the highest reflectance (percentile = 100).  If this
% value is < 1, we do nothing.  If the value exceeds one, we scale the
% illuminant so that the max reflectance is 1.
%
% If there are specularities and this scaling makes the reflectances too
% small.  To allow for these cases, choose a percentile less than 100, say
% 98.
%
% Example:
%  fullFileName = fullfile(isetRootPath,'data','images', ...
%                   'multispectral','StuffedAnimals_tungsten-hdrs.mat');
%  scene = sceneFromFile(fullFileName,'multispectral');
%  scene = sceneAdjustReflectance(scene);
%  vcAddAndSelectObject('scene',scene); sceneWindow;
%
%  scene = sceneAdjustReflectance(scene,98);
%  vcAddAndSelectObject('scene',scene); sceneWindow;
%
% See also:  sceneAdjustIlluminant, sceneAdjustLuminance
%
% (c) Imageval, 2012

if notDefined('scene'), error('Scene required.'); end
if notDefined('prct'),  prct = 100; end

% Find the current peak reflectance value
r = sceneGet(scene, 'reflectance');

% Find the value that we want to be at unit reflectance.  Because of
% specularities, it may be that some locations will have > 1 reflectance.
if prct == 100, peakReflectance = max(r(:));
else
    % Faster band-wise. Not quite the same as the global percentile. This
    % is the band-wise percentile.
    nWave = sceneGet(scene,'nwave');
    peakReflectanceBand = zeros(nWave,1);
    str = sprintf('Finding %d percentile reflectance',prct);
    showBar = ieSessionGet('waitbar');
    if showBar, h = waitbar(0,str); end
    for ii=1:nWave
        if showBar, waitbar(ii/nWave,h); end
        tmp = r(:,:,ii);
        peakReflectanceBand(ii) = iePrctile(tmp(:),prct);
    end
    if showBar, close(h); end
    peakReflectance = iePrctile(peakReflectanceBand,prct);
end

if peakReflectance > 1
    newIlluminant = sceneGet(scene,'illuminant photons')*peakReflectance;
    scene = sceneSet(scene,'illuminant photons',newIlluminant);
else
    disp('All reflectance values are < 1');
end

end