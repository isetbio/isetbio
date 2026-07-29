%% Test roundtrip spectra

%% Initialize
ieInit;
if ~mcDockerExists, mcDockerConfig; end % check whether we can use docker

%% Generate correct xyz2rgb matrix
d = displayCreate(); % We (can) use these primaries in PBRT
rgb2xyz = displayGet(d,'rgb2xyz');
xyz2rgb = inv(rgb2xyz);

wave = displayGet(d,'wave');
rgbPrimaries = displayGet(d,'spd primaries');
maxB = max(rgbPrimaries(:,3));
load('RGBPrimaries2.mat');
rgbPrimaries = [R G B];
% rgbPrimaries = rgbPrimaries.*(maxB/max(rgbPrimaries(:,3)));
figure(); hold on; grid on;
plot(wave,rgbPrimaries(:,1),'r');
plot(wave,rgbPrimaries(:,2),'g');
plot(wave,rgbPrimaries(:,3),'b');

%% Test something
%{
expectedEnergy = [1 1 1] * rgbPrimaries';
expectedPhotons = Energy2Quanta(wave,expectedEnergy');
figure();
plot(wave,expectedPhotons,'b'); hold on;

load('r.mat');
% r = r.*(max(expectedPhotons)/max(r));
plot(wave,r,'r');

expectedXYZ = ieXYZFromEnergy(expectedEnergy,wave);
expectedXYZ = expectedXYZ*(100/expectedXYZ(2));
fprintf('ExpectedXYZ: %0.2f %0.2f %0.2f \n',expectedXYZ)
expectedRGB = imageLinearTransform(expectedXYZ, xyz2rgb);
fprintf('ExpectedRGB: %0.2f %0.2f %0.2f \n',expectedRGB)
%}

%% Render a few RGB pixels of known value

im = ones(10,10,3);
imshow(im); title('Original Scene');

imDir = fullfile(isetbioRootPath,'local','imageTextures');
if(~exist(imDir))
    mkdir(imDir);
end

imFile = fullfile(imDir,'GT_image.png');
imwrite(im,imFile);

RGB1 = im2double(imread(imFile));

distance = 1;
planeFOV = 20;
width = 2*tand(planeFOV/2)*distance;
sz = [width width];

% Try with ISET
%{
scene = sceneFromFile(imFile,'rgb');
% scene = sceneAdjustIlluminant(scene,'equalEnergy.mat');
scene = sceneSet(scene,'distance',distance);
scene = sceneSet(scene,'wangular',planeFOV);
meanLum = sceneGet(scene,'mean luminance');
ieAddObject(scene);
sceneWindow;
xyz = sceneGet(scene,'xyz');
%}
%{
photons = piReadDAT('/Users/tlian/GitRepos/pbrt-v3-spectral/Debug/pbrt.dat');
ieObj1 = piOICreate(photons);

scenePBRT = scene;
scenePBRT = sceneAdjustIlluminant(scenePBRT,'equalEnergy.mat');
scenePBRT = sceneSet(scenePBRT,'wave',sceneGet(oi,'wave'));
scenePBRT = sceneSet(scenePBRT,'photons',photons);
scenePBRT = sceneSet(scenePBRT,'mean luminance',meanLum);
ieAddObject(scenePBRT);
sceneWindow;
%}

se1 = sceneEye('texturedPlane',...
    'planeDistance',distance,...
    'planeSize',sz,...
    'planeTexture',imFile,...
    'useDisplaySPD',1);

se1.fov = planeFOV;

se1.resolution = 10;
se1.numRays = 64;
se1.numBounces = 1;
se1.numCABands = 0;
se1.accommodation = 1;

se1.debugMode = true;
se1.name = 'first_pass';

[ieObj1,~,sf] = se1.render();

% Convert units
wave = sceneGet(ieObj1,'wave');
%{
energy = sceneGet(ieObj1,'photons');
photons = Energy2Quanta(wave,energy);
ieObj1 = sceneSet(ieObj1,'photons',photons);
ieObj1 = sceneAdjustLuminance(ieObj1, 100);
%}

ieAddObject(ieObj1);
sceneWindow;

photons = sceneGet(ieObj1,'photons');
pbrtPhotons = squeeze(mean(mean(photons,1),2));

spectraFigure = figure(); hold on;
plot(wave,pbrtPhotons);
%{
% Plot range of spectra
sz = sceneGet(ieObj1,'size');
for ii = 1:sz(1)
    for jj = 1:sz(2)
        plot(wave,squeeze(photons(ii,jj,:)));
    end
end
%}
grid on;
xlabel('Wavelength (nm)');
ylabel('Radiance')

%{
% Check SPD
% Technically the output from PBRT is energy since we are using the RGB
% primaries.
% expectedEnergy = [1 1 1] * rgbPrimaries'*sf*6.443*0.0036649; % Still not sure where the 6.443 comes from.
expectedEnergy = [1 1 1] * rgbPrimaries'*0.0036649*6.443; 
% expectedPhotons = expectedPhotons.*(max(pbrtPhotons)/max(expectedPhotons));
expectedPhotons = Energy2Quanta(expectedEnergy,wave);
plot(wave,expectedPhotons,'r--');

expectedXYZ = ieXYZFromPhotons(expectedPhotons',wave);
expectedXYZ = expectedXYZ*(100/expectedXYZ(2));
fprintf('ExpectedXYZ: %0.2f %0.2f %0.2f \n',expectedXYZ)
expectedRGB = imageLinearTransform(expectedXYZ, xyz2rgb);
fprintf('ExpectedRGB: %0.2f %0.2f %0.2f \n',expectedRGB)

XYZ2 = sceneGet(ieObj1,'xyz');
measuredXYZ = squeeze(mean(mean(XYZ2,1),2))';
fprintf('MeasuredXYZ: %0.2f %0.2f %0.2f \n',measuredXYZ)
measuredRGB = imageLinearTransform(measuredXYZ, xyz2rgb);
fprintf('MeasuredRGB: %0.2f %0.2f %0.2f \n',measuredRGB)
%}

% Check range of XYZ/RGB values
%{
for ii = 1:sz(1)
    for jj = 1:sz(2)
        measuredXYZ = squeeze(XYZ2(ii,jj,:))';
        measuredRGB = imageLinearTransform(measuredXYZ, xyz2rgb);
        fprintf('%0.2f %0.2f %0.2f \n',measuredRGB)
    end
end
%}

% Check deltaE (needs isetCam)
%{
whitePnt = displayGet(d,'whitepoint');
dEab = deltaEab(measuredXYZ,expectedXYZ,whitePnt);
fprintf('dEab: %0.2f \n',dEab)
%}

%% Save out a new RGB file, repeat.

XYZ2 = sceneGet(ieObj1,'xyz');
RGB2 = imageLinearTransform(XYZ2, xyz2rgb);

imFile = fullfile(isetbioRootPath,'local','imageTextures','display_image.png');
imwrite(RGB2,imFile);

distance = 1;
planeFOV = 20;
width = 2*tand(planeFOV/2)*distance;
sz = [width width];
    
se2 = sceneEye('texturedPlane',...
    'planeDistance',distance,...
    'planeSize',sz,...
    'planeTexture',imFile,...
    'useDisplaySPD',1);

se2.fov = planeFOV;

se2.resolution = 10;  
se2.numRays = 64;
se2.numBounces = 1;
se2.numCABands = 0;
se2.name = 'second_pass';

se2.accommodation = 1;

se2.debugMode = true;

ieObj2 = se2.render();

wave = sceneGet(ieObj2,'wave');
%{
energy = sceneGet(ieObj2,'photons');
photons = Energy2Quanta(wave,energy);
ieObj2 = sceneSet(ieObj2,'photons',photons);
ieObj2 = sceneAdjustLuminance(ieObj2, 100);
%}

ieAddObject(ieObj2);
sceneWindow;

XYZ3 = sceneGet(ieObj2,'xyz');
RGB3 = imageLinearTransform(XYZ3, xyz2rgb);


fprintf('RGB 1: [%0.2f %0.2f %0.2f]\n',squeeze(mean(mean(RGB1,1),2)))
fprintf('RGB 2: [%0.2f %0.2f %0.2f]\n',squeeze(mean(mean(RGB2,1),2)))
fprintf('RGB 3: [%0.2f %0.2f %0.2f]\n',squeeze(mean(mean(RGB3,1),2)))

photons = sceneGet(ieObj2,'photons');
wave = sceneGet(ieObj2,'wave');
pbrtPhotons = squeeze(mean(mean(photons,1),2));

figure(spectraFigure);
plot(wave,pbrtPhotons);

legend('Pass 1','Pass 2');

