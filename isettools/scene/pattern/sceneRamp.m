function scene = sceneRamp(scene,sz)
%% Intensity ramp (see L-star chart for L* steps)

if notDefined('sz'), sz = 128; end

scene = sceneSet(scene,'name','ramp');
scene = initDefaultSpectrum(scene,'hyperspectral');
nWave = sceneGet(scene,'nwave');
wave = sceneGet(scene,'wave');

img = imgRamp(sz);
img = img/(max(img(:)));

il = illuminantCreate('equal photons',wave,100);
scene = sceneSet(scene,'illuminant',il);

img = repmat(img,[1,1,nWave]);
[img,r,c] = RGB2XWFormat(img);
illP = illuminantGet(il,'photons');
img = img*diag(illP);
img = XW2RGBFormat(img,r,c);
scene = sceneSet(scene,'photons',img);

end