%% s_sceneIlluminantSpace
%
% Illustrate the spatial spectral illumination representation
%
% (c) Imageval Consulting, LLC 2012

%%
s_initISET

%%
scene = sceneCreate;

%% Make a version of the illuminant that matches the size of the scene.

illE = sceneGet(scene,'illuminant energy');
sz = sceneGet(scene,'size');
r = sz(1); c = sz(2);
illE2 = repmat(illE(:),1,r*c);
illE2 = XW2RGBFormat(illE2',r,c);
scene = sceneSet(scene,'illuminant energy',illE2);

%% Have a look at the scene illuminant - it will look blank/white
%
illEnergy = sceneGet(scene,'illuminant energy');
[illEnergy,r,c] = RGB2XWFormat(illEnergy);
wave = sceneGet(scene,'wave');
XYZ = ieXYZFromEnergy(illEnergy,wave);
srgb = xyz2srgb(XW2RGBFormat(XYZ,r,c));
vcNewGraphWin; imagesc(srgb);  % Looks white, I guess

%% Make a spatial spectral illuminant matrix  
% We adjust the SPD along the rows to go from 6500 to 3000 K
illPhotons = sceneGet(scene,'illuminant photons');
[r,c,w] = size(illPhotons);

cTemp = linspace(6500,3000,r);
spd   = blackbody(wave,cTemp);
% vcNewGraphWin; plot(wave,spd);

% START HERE
for rr=1:r
    illPhotons(rr,:,:) = squeeze(illPhotons(rr,:,:)) * diag((spd(:,rr)./illE(:)));
end

% Correct the energy for this illuminant
reflectance = sceneGet(scene,'reflectance');
p = reflectance .* illPhotons;

% When we divide to obtain the reflectance it should be the same
scene = sceneSet(scene,'photons',p);
scene = sceneSet(scene,'illuminant photons',illPhotons);

vcAddAndSelectObject(scene); sceneWindow;

%% Show the result as a graph
illEnergy = sceneGet(scene,'illuminant energy');
vcNewGraphWin;
mesh(wave,1:r,squeeze(illEnergy(:,1,:)))
xlabel('wavelength')
ylabel('pos')

%% Show the illuminant color change as an image
illEnergy = sceneGet(scene,'illuminant energy');
[illEnergy,r,c] = RGB2XWFormat(illEnergy);
wave = sceneGet(scene,'wave');
XYZ = ieXYZFromEnergy(illEnergy,wave);
XYZ = XW2RGBFormat(XYZ,r,c);
srgb = xyz2srgb(XYZ);
vcNewGraphWin; imagesc(srgb);  % Looks white, I guess

%% Original scene for comparison
sceneM = sceneCreate;
vcAddAndSelectObject(sceneM); sceneWindow;

%% Check conversion
sceneGet(scene,'illuminant format')

illPhotons = sceneGet(scene,'illuminant photons');
scene = sceneSet(scene,'illuminant photons',illPhotons);
vcAddAndSelectObject(scene); sceneWindow;

%% Make an intensity varying illuminant, starting with spatial spectral
% We adjust the intensity along the cols 
illPhotons = sceneGet(scene,'illuminant photons');
[r,c,w] = size(illPhotons);

cc = 1:c;
illScale = 1 + 0.5*sin(2*pi*(cc/c));
% vcNewGraphWin; plot(wave,spd);

for cc=1:c
    illPhotons(:,cc,:) = squeeze(illPhotons(:,cc,:)) * illScale(cc);
end

% Correct the energy for this illuminant
reflectance = sceneGet(scene,'reflectance');
p = reflectance .* illPhotons;

% When we divide to obtain the reflectance it should be the same
scene = sceneSet(scene,'photons',p);
scene = sceneSet(scene,'illuminant photons',illPhotons);

vcAddAndSelectObject(scene); sceneWindow;

%% Now do the same, but starting with a spectral and scaling across cols
scene = sceneCreate;
vcAddAndSelectObject(scene); sceneWindow;

illP = sceneGet(scene,'illuminant photons');
sz = sceneGet(scene,'size');
r = sz(1); c = sz(2);
illP2 = repmat(illP(:),1,r*c)';
illP2 = XW2RGBFormat(illP2,r,c);

cc = 1:c;
freq = 1;
illScale = 1 + 0.5*sin(2*pi*freq*(cc/c));
vcNewGraphWin; plot(illScale); grid on

% Scale across the columns
for cc=1:c
    illP2(:,cc,:) = squeeze(illP2(:,cc,:)) * illScale(cc);
end

% Correct the energy for this illuminant
reflectance = sceneGet(scene,'reflectance');
scene = sceneSet(scene,'illuminant photons', illP2);
scene = sceneSet(scene,'photons',illP2 .* reflectance);

vcAddAndSelectObject(scene); sceneWindow;

%%
