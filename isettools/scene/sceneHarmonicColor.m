function [scene,p] = sceneHarmonicColor(scene,params, wave)
%% Create a scene of a (windowed) harmonic function with RGB components.
%
% Harmonic parameters are: parms.freq, parms.row, parms.col, parms.ang
% parms.ph, parms.contrast
%
% Missing default parameters are supplied by imageHarmonic.
%
% The frequency is with respect to the image (cyces/image).  To determine
% cycles/deg, use cpd: freq/sceneGet(scene,'fov');
%

scene = sceneSet(scene,'name','harmonic');

if notDefined('wave')
    scene = initDefaultSpectrum(scene,'hyperspectral');
else
    scene = initDefaultSpectrum(scene, 'custom',wave);
end

nWave = sceneGet(scene,'nwave');

% TODO: Adjust pass the parameters back from the imgHarmonic window. In
% other cases, they are simply attached to the global parameters in
% vcSESSION.  We can get them by a getappdata call in here, but not if we
% close the window as part of imageSetHarmonic
if notDefined('params')
    [h, params] = imageSetHarmonic; waitfor(h);
    img = imageHarmonicColor(params);
    p   = params;
else
    [img,p] = imageHarmonicColor(params);
end


scene = sceneFromFile(stimulusRGBdata, 'rgb', params.meanLuminance, display);

% To reduce rounding error problems for large dynamic range, we set the
% lowest value to something slightly more than zero.  This is due to the
% ieCompressData scheme.
img(img==0) = 1e-4;
img   = img/(2*max(img(:)));    % Forces mean reflectance to 25% gray

% Mean illuminant at 100 cd
wave = sceneGet(scene,'wave');
il = illuminantCreate('equal photons',wave,100);
scene = sceneSet(scene,'illuminant',il);

img = repmat(img,[1,1,nWave]);
[img,r,c] = RGB2XWFormat(img);
illP = illuminantGet(il,'photons');
img = img*diag(illP);
img = XW2RGBFormat(img,r,c);
scene = sceneSet(scene,'photons',img);

% set scene field of view
scene = sceneSet(scene, 'h fov', 1);

end
