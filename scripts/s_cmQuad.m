%% Show the basic quadrature calculation
%
% This is for a 1D image.
% Tomorrow extend the calculations to a 2D image.
%

%%
ieInit

%% Testing the basic idea

%  Make a signal and two harmonics in quadrature
x = [0:127]/128;
f = 2;
s = sin(2*pi*f*x);
c = cos(2*pi*f*x);
sig = 0.4*square(2*f*pi*x) + 0.5;

eBase = dot(sig,s)^2 + dot(sig,c)^2;

%
vcNewGraphWin; plot(x,sig,'k-',x,s,'r-',x,c,'b-');

%% Take the inner product of the signal with each harmonic. 

% Then compute the energy, also known as the amplitude at that
% frequency.

% Shift the signal and recompute
for ii=1:2:10
    eShift = dot(circshift(sig,ii),s)^2 + dot(circshift(sig,ii),c)^2;
    fprintf('Difference: %.6f\n',eBase - eShift)
end

%% Now do the same, but for a 2D image
img    = repmat(sig,[128,1]);
simg   = repmat(s,[128,1]);
cimg   = repmat(c,[128,1]);
vcNewGraphWin; imagesc(img); colormap(gray); axis image

eBase = dot(img(:),simg(:))^2 + dot(img(:),cimg(:))^2;
for ii=1:2:10
    thisIMG = circshift(img,ii,2);
    imagesc(thisIMG); colormap(gray); axis image; pause(0.2);
    eShift = dot(thisIMG(:),simg(:))^2 + dot(thisIMG(:),cimg(:))^2;
    fprintf('Difference: %.6f\n',eBase - eShift)
end

%% Now modify the calculation by applying by a Gaussian envelope

% Big difference with a small envelope, and little difference with a big
% envelope, like the full harmonic above.
spread = 32;
g = fspecial('gaussian',[128 128],spread);

vcNewGraphWin; imagesc(g); colormap(gray); axis image
gsimg = g .* simg;
gcimg = g .* cimg;

eBase = dot(img(:),gsimg(:))^2 + dot(img(:),gcimg(:))^2;

% Small shifts and there is a constant response
fprintf('\n----\n');
for ii=1:2:spread
    thisIMG = circshift(img,ii,2);
    eShift = dot(thisIMG(:),gsimg(:))^2 + dot(thisIMG(:),gcimg(:))^2;
    fprintf('Difference (percentage): %.6f (step %d)\n',100*(eBase - eShift)/eBase,ii)
end

%%  Now try with some stimuil

ieInit
chdir(fullfile(isetbioRootPath,'local'));

%% Harmonic calculations

% In this case, the pattern is genuinely 1D and although the eye movements
% are in both directions, the filters are 1D, and the noise will be varied.

hparams = harmonicP; 
hparams.freq = 3;
scene = sceneCreate('harmonic',hparams);
fov = 4;
scene = sceneSet(scene,'fov',3);
oi = oiCreate; oi = oiCompute(oi,scene);
ieAddObject(oi);

%% Start the cone absorptions with no Poisson noise

cm = coneMosaic;
cm.setSizeToFOV(fov*0.7);
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.noiseFlag = 'none';
cm.compute(oi);
[r,c,t] = size(cm.absorptions);
% cm.window;

%% Calculate for the quadrature and a true circular shift

% Make the Gabor patches in quadrature phase as above.
hparams = harmonicP;
hparams.row = size(cm.absorptions,1);
hparams.col = size(cm.absorptions,2);
hparams.freq = 2;
hparams.ph = pi/2;
sQuad = imageHarmonic(hparams);
sQuad = sQuad - 1;
% vcNewGraphWin; mesh(sQuad); colormap(gray)

hparams.ph = 0;
cQuad = imageHarmonic(hparams);
cQuad = cQuad - 1;
% vcNewGraphWin; mesh(cQuad); colormap(gray)

%% Use this as the base image and just shift it

baseFrame   = 2;
baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% Try just shifting the base frame.

% We get a very solid result.  No noise and pure shifting has no impact
for ii=1:t
    thisIMG = circshift(baseIMG,ii,2);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
%{
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');
%}

%% Instead of shifting, use the eye movement sequence.

% In this case, the cone images are not shifted copies of one another.  And
% when the eye movement shift is large the change is rather significant.
for ii=1:t
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');

%%  Does the failure arise in part because of the different types of cones

% So, here is a mosaic with only M cones.   Still no noise
cm.spatialDensity = [0 0 1 0];
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
[r,c,t] = size(cm.absorptions);
% cm.window;
% cm.plot('eye movement path')
baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% First try just shifting the base frame.  In the no noise case, we get a
% very solid result.
for ii=1:t
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');

%% Add in the Poisson noise

cm.spatialDensity = [0 0 1 0];
cm.noiseFlag = 'random';
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
[r,c,t] = size(cm.absorptions);
% cm.window;

baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% First try just shifting the base frame.  In the no noise case, we get a
% very solid result.
for ii=1:t
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');

%% And now half L and half M cones

cm.spatialDensity = [0 .5 0.5 0];
cm.noiseFlag = 'random';
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
[r,c,t] = size(cm.absorptions);
% cm.window;

baseIMG = cm.absorptions(:,:,baseFrame);
eBase   = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% First try just shifting the base frame.  In the no noise case, we get a
% very solid result.
for ii=1:t
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');

%% Finally, put the S cones and noise back in

cm.spatialDensity = [0 .6 0.3 0.1];
cm.noiseFlag = 'random';
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
[r,c,t] = size(cm.absorptions);
% cm.window;

%
baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% First try just shifting the base frame.  In the no noise case, we get a
% very solid result.
for ii=1:t
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');
%% ----------- Dead leaves calculation ----------------

% In this case, the pattern is 2D. 
% The filters are still 1D

dlSize = 128; dlSigma = 3;
n = 128; sigma = 3; 
scene = sceneDeadleaves(dlSize, dlSigma);
fov = 4;
scene = sceneSet(scene,'fov',3);
% sceneWindow(scene);

oi = oiCreate; oi = oiCompute(oi,scene);
ieAddObject(oi);
% oiWindow(oi);

%% Start the cone absorptions with no Poisson noise

cm = coneMosaic;
cm.setSizeToFOV(fov*0.7);
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.noiseFlag = 'random';
cm.compute(oi);
[r,c,t] = size(cm.absorptions);
cm.plot('eye movement path');
% cm.window;

%% Calculate for the quadrature and a true circular shift

% Make the Gabor patches in quadrature phase as above.
hparams = harmonicP;
hparams.row = size(cm.absorptions,1);
hparams.col = size(cm.absorptions,2);
hparams.freq = 4;
hparams.ph = pi/2;
sQuad = imageHarmonic(hparams);
sQuad = sQuad - 1;
% vcNewGraphWin; mesh(sQuad); colormap(gray)

hparams.ph = 0;
cQuad = imageHarmonic(hparams);
cQuad = cQuad - 1;
% vcNewGraphWin; mesh(cQuad); colormap(gray)

%% Use this as the base image and just shift it

baseFrame   = 2;
baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

%{
% Try just shifting the base frame.
% We get a very solid result.  No noise and pure shifting has no impact
% This was so repeatable, I commented it out.
for ii=1:t
    thisIMG = circshift(baseIMG,ii,2);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end
%}

%{
% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');
%}

%% Now use the absorptions

for ii=1:t
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end
% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');

%% Compute the photocurrent

cm.computeCurrent;
baseFrame   = 25;
baseIMG     = cm.current(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

%% Compute the error, but because of the ramp up use 25 as the base frame
for ii=1:t
    thisIMG     = cm.current(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end
% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');

%%
