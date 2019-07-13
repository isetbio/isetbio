%% Show the quadrature calculation at the cone excitations for a harmonic
%
% On the cone excitations. When the shift is within a few cones the
% response stays within a couple of percent. But for a shift of five
% or so cones, the value changes more than a couple of percent.
%

%% Quadrature filters with cone mosaic data

ieInit
chdir(fullfile(isetbioRootPath,'local'));

%% Harmonic image at the cone mosaic

% In this case, the pattern is 1D and although the eye movements are
% in both directions, the filters are 1D, and the noise will be
% varied.
hparams = harmonicP; 
hparams.freq = 3;
scene = sceneCreate('harmonic',hparams);
fov = 4; scene = sceneSet(scene,'fov',3);
oi = oiCreate; oi = oiCompute(oi,scene); ieAddObject(oi);

cm = coneMosaic;
cm.noiseFlag = 'none';
cm.spatialDensity = [0 .6 .3 .1];
cm.setSizeToFOV(fov*0.7);
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
% cm.window;

%% Calculate the quadrature filters

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

%% Validate that the circshift does the expected 0 change

baseFrame   = 2;
baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(cm.tSamples,1);

% Try just shifting the base frame. We get a clean result.  No noise
% and pure shifting has no impact
for ii=1:cm.tSamples
    thisIMG = circshift(baseIMG,ii,2);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
% {
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on; set(gca, 'ylim',[-1 1]); title('Circular shift')
xlabel('Distance (cones)'); ylabel('Percent error');
%}

%% Instead of shifting, use the eye movement sequence.

% There is still no Poisson noise.  So frame-to-frame is just a shift
% but no independent noise samples.

% The cone images are not shifted copies of one another because there
% is some clipping at the edges;  new signals moving in and out. When
% the eye movement shift is large the change is rather significant.
%
% For a 64 x 64 image, a 10 cone horizontal shift changes 20 percent
% of the columns.
for ii=1:cm.tSamples
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
subplot(1,2,1), plot(d,eAbsorb,'o-'); grid on; set(gca,'ylim',[-10 10]);
xlabel('Distance (cones)'); ylabel('Percent error'); 
subplot(1,2,2), cm.plot('eye movement path','hf',gca);
title(sprintf('Noise %s',cm.noiseFlag));

%% Recompute the cone mosaic, adding Poisson noise

cm.spatialDensity = [0 .6 .3 .1];
cm.noiseFlag = 'random';
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
% cm.window;

baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% First try just shifting the base frame.  In the no noise case, we get a
% very solid result.
for ii=1:cm.tSamples
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
subplot(1,2,1), plot(d,eAbsorb,'o-'); grid on; set(gca,'ylim',[-10 10]);
xlabel('Distance (cones)'); ylabel('Percent error');
subplot(1,2,2), cm.plot('eye movement path','hf',gca);
title(sprintf('Noise %s',cm.noiseFlag));

%% END

%{
%% Same pattern with only one cone type

cm.spatialDensity = [0 0 1 0];
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
% cm.window;

baseIMG = cm.absorptions(:,:,baseFrame);
eBase   = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% First try just shifting the base frame.  In the no noise case, we get a
% very solid result.
for ii=1:cm.tSamples
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
subplot(1,2,1), plot(d,eAbsorb,'o-'); grid on; set(gca,'ylim',[-10 10]);
xlabel('Distance (cones)'); ylabel('Percent error');
subplot(1,2,2), cm.plot('eye movement path','hf',gca);
title(sprintf('Noise %s',cm.noiseFlag));
%}


%{
%%  Does the failure arise in part because of the different types of cones

% So, here is a mosaic with only M cones.   Still no noise
cm.spatialDensity = [0 0 1 0];
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);

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
subplot(1,2,1), plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');
subplot(1,2,2), cm.plot('eye movement path','hf',gca);

%% Finally, put the S cones and noise back in

cm.spatialDensity = [0 .6 0.3 0.1];
cm.noiseFlag = 'random';
cm.emGenSequence(50,'rSeed',[],'nTrials',1);
cm.compute(oi);
% cm.window;

%
baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(t,1);

% First try just shifting the base frame.  In the no noise case, we get a
% very solid result.
for ii=1:cm.tSamples
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
subplot(1,2,1), plot(d,eAbsorb,'o-'); grid on;
xlabel('Distance (cones)'); ylabel('Percent error');
subplot(1,2,2), cm.plot('eye movement path','hf',gca);


%}
