function [scene,params] = sceneHarmonic(scene, params, wave)
%% Create a scene of a (windowed) harmonic function.
%
% Syntax
%  [scene,p] = sceneHarmonic(scene,params, wave)
%
% Description
%  Create a Gabor scene (harmonic modulated by a Gaussian). The Harmonic
%  parameters are:  
%   params.name 
%   params.ang
%   params.contrast
%   params.freq
%   params.ph
%   params.row, parms.col
%   params.GaborFlag
%
%  If you do not supply harmonic parameters, we use the default returned by
%  the function harmonicP
%
% The frequency is with respect to the image (cyces/image).  To determine
% cycles/deg, use cpd: freq/sceneGet(scene,'fov');
%
% Imageval Consulting, LLC Copyright 2006

% Examples
%{
  params = harmonicP;
  scene = sceneInit;
  scene = sceneHarmonic(scene,params);
  scene = sceneComplete(scene);
  ieAddObject(scene); sceneWindow;
%}
%{
  params.freq = 5;
  params.ang  = pi/3;
  params.GaborFlag = 0.2;
  scene = sceneHarmonic(scene,params);
  ieAddObject(scene); sceneWindow;
%}

%% Build the image from harmonic parameters
if notDefined('params')
    params = harmonicP; 
    warning('Using default harmonic parameters'); 
end
img = imageHarmonic(params);

% To reduce rounding error problems for large dynamic range, we set the
% lowest value to something slightly more than zero.
img(img==0) = 1e-4;
img   = img/(2*max(img(:)));    % Forces mean reflectance to 25% gray

%% Build scene from img 

if notDefined('scene'), scene = sceneInit; end
scene = sceneSet(scene,'name','harmonic');
if notDefined('wave')
    scene = initDefaultSpectrum(scene,'hyperspectral');
else
    scene = initDefaultSpectrum(scene, 'custom',wave);
end

nWave = sceneGet(scene,'nwave');

% Mean illuminant at 100 cd
wave  = sceneGet(scene,'wave');
il    = illuminantCreate('equal photons',wave,100);
scene = sceneSet(scene,'illuminant',il);

% Build up the photons
photons = repmat(img,[1,1,nWave]);
[photons,r,c] = RGB2XWFormat(photons);
illP = illuminantGet(il,'photons');
photons = photons*diag(illP);
photons = XW2RGBFormat(photons,r,c);
scene = sceneSet(scene,'photons',photons);

%% set scene field of view
scene = sceneSet(scene, 'h fov', 1);
scene = sceneAdjustLuminance(scene,100);

end