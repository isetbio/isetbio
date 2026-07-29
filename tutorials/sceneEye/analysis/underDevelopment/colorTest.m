%% Test roundtrip spectra

%% Initialize
ieInit;
if ~mcDockerExists, mcDockerConfig; end % check whether we can use docker

load('pbrtOutputUpdated2.mat');
%% Generate correct xyz2rgb matrix

d = displayCreate(); % We use these primaries in PBRT
d = displaySet(d,'wave',400:10:700);

rgb2xyz = displayGet(d,'rgb2xyz');
xyz2rgb = inv(rgb2xyz).*(1/0.0088);

wave = displayGet(d,'wave');
rgbPrimaries = displayGet(d,'spd primaries');
maxB = max(rgbPrimaries(:,3));
rgbPrimaries = rgbPrimaries./maxB;

figure(); hold on; grid on;
plot(wave,rgbPrimaries(:,1),'r');
plot(wave,rgbPrimaries(:,2),'g');
plot(wave,rgbPrimaries(:,3),'b');

plot(wave,R,'r:');
plot(wave,G,'g:');
plot(wave,B,'b:');

% Test case
testRGB = [0 0 1]; % gamma corrected
testRGB_linear = testRGB.^2.2;
spectrumEnergy = testRGB_linear*rgbPrimaries'; % RGB to spectrum (energy)
% plot(wave,spectrumEnergy);
xyz = ieXYZFromEnergy(spectrumEnergy,wave); % spectrum to XYZ
xyz = xyz/max(xyz(:)); % Normalize
calcRGB_linear = imageLinearTransform(reshape(xyz,[1 1 3]),xyz2rgb.*(1/0.0088));
calcRGB = calcRGB_linear.^(1/2.2);

fprintf('testRGB = [%f %f %f]\n',testRGB(1),testRGB(2),testRGB(3));
fprintf('calcRGB = [%f %f %f]\n',calcRGB(1),calcRGB(2),calcRGB(3));

%% Render a few RGB pixels of known value
imSize = 512;

% A blue patch
rgbOrig = [0 0 1];
imOrig = zeros(imSize,imSize,3);
imOrig(:,:,3) = 1; % We assume this is lrgb

% imshow(imOrig); title('Original RGB');

imDir = fullfile(isetbioRootPath,'local','imageTextures');
if(~exist(imDir))
    mkdir(imDir);
end

imFile = fullfile(imDir,'GT_image.png');
imwrite(imOrig,imFile);

RGB1 = im2double(imread(imFile));

distance = 1;
planeFOV = 20;
width = 2*tand(planeFOV/2)*distance;
sz = [width width];

se1 = sceneEye('texturedPlane',...
    'planeDistance',distance,...
    'planeSize',sz,...
    'planeTexture',imFile,...
    'useDisplaySPD',1,...
    'gamma','false');

se1.fov = planeFOV;

se1.resolution = imSize;
se1.numRays = 64;
se1.numBounces = 1;
se1.numCABands = 0;
se1.accommodation = 1;

se1.debugMode = true;
se1.name = 'renderedRGB';

[ieObj1,~,sf] = se1.render();

%% Compare photons

load('pbrtOutputUpdated2.mat');
energyDirectPBRT = L;

photonsRendered = sceneGet(ieObj1,'photons');
sz = size(photonsRendered);
midPt = floor(sz(1)/2);
photonsRendered = squeeze(photonsRendered(midPt,midPt,:));
energyRendered = Quanta2Energy(wave, photonsRendered);
%{
wave = sceneGet(ieObj1,'wave');
scaleFactor = energyDirectPBRT(wave == 450)/...
    energyRendered(wave == 450);
energyRendered = energyRendered.*scaleFactor;

figure(); hold on;
plot(wave,energyRendered);
plot(wave,energyDirectPBRT);
plot(wave,spectrumEnergy);
xlim([min(wave) max(wave)])
legend('Rendered','Rendered(direct)','Calculated');
%}

%% Compare conversion to RGB

% rendered
xyz = ieXYZFromEnergy(energyRendered,wave); % spectrum to XYZ
xyz = xyz/max(xyz(:)); % Normalize
renderedlRGB = imageLinearTransform(reshape(xyz,[1 1 3]),xyz2rgb);
renderedlRGB = max(renderedlRGB,0);
renderedRGB = renderedlRGB.^(1/2.2);

% calculated

xyz = ieXYZFromEnergy(spectrumEnergy,wave); % spectrum to XYZ
xyz = xyz/max(xyz(:)); % Normalize
calclRGB = imageLinearTransform(reshape(xyz,[1 1 3]),xyz2rgb);
calcRGB = calclRGB.^(1/2.2);

fprintf('renderedRGB = [%f %f %f]\n',renderedRGB(1),renderedRGB(2),renderedRGB(3));
fprintf('calcRGB = [%f %f %f]\n',calcRGB(1),calcRGB(2),calcRGB(3));


%% Compare final RGB values

% Crop
cropPixels = 3;
ieObj1 = oiCropBorder(ieObj1,cropPixels);
imOrig =  imOrig(cropPixels:size(imOrig,1) - cropPixels,...
    cropPixels:size(imOrig,2) - cropPixels,:);

% Get RGB
imRendered = getDisplayRGB(ieObj1);
% imRendered = imRendered.*(2.2);
rgbRendered = squeeze(imRendered(midPt,midPt,:)); % center pixel

figure();
subplot(1,2,1); imshow(imOrig);
title('Original')
subplot(1,2,2); imshow(imRendered);
title('Rendered')

% Compare RGB
fprintf('rgbOrig = [%f %f %f]\n',rgbOrig(1),rgbOrig(2),rgbOrig(3));
fprintf('rgbRendered = [%f %f %f]\n',rgbRendered(1),rgbRendered(2),rgbRendered(3));

%% MSSSIM
pxPerDeg = 25.6;
[weights, freq] = getMSSSIMweights(100,pxPerDeg,0);
figure();
for ii = 1:3
[msssimvals(ii), msssimmap] = msssim(im2uint8(imOrig(:,:,ii)),...
    im2uint8(imRendered(:,:,ii)),weights);
subplot(1,3,ii);
imagesc(msssimmap); colormap(gray); colorbar;
axis image; axis off;

end

