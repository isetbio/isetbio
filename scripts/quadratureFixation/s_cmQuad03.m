%% ----------- Dead leaves calculation ----------------
%
%
%
%%
ieInit
chdir(fullfile(isetbioRootPath,'local'));

%% Make the dead leaves scene, oi and cone responses with Poisson noise

dlSize = 128; dlSigma = 3;
n = 128; sigma = 3; 
scene = sceneDeadleaves(dlSize, dlSigma);
fov = 4; scene = sceneSet(scene,'fov',3);
% sceneWindow(scene);

oi = oiCreate; oi = oiCompute(oi,scene); ieAddObject(oi);
% oiWindow(oi);

cm = coneMosaic;
cm.setSizeToFOV(fov*0.7);
cm.emGenSequence(100,'rSeed',[],'nTrials',1);
cm.noiseFlag = 'random';
cm.compute(oi);
% cm.plot('eye movement path');
% cm.window;

%% Calculate the quadrature filters

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
eAbsorb = zeros(cm.tSamples,1);

% Try just shifting the base frame.
% We get a very solid result.  No noise and pure shifting has no impact
% This was so repeatable, I commented it out.
for ii=1:t
    thisIMG = circshift(baseIMG,ii,2);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% {
% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on; set(gca, 'ylim',[-1 1]); title('Circular shift')
xlabel('Distance (cones)'); ylabel('Percent error');
%}

%% Now use the absorptions

for ii=1:cm.tSamples
    thisIMG     = cm.absorptions(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
subplot(1,2,1), plot(d,eAbsorb,'o-'); grid on; set(gca,'ylim',[-10 10]);
xlabel('Distance (cones)'); ylabel('Percent error');
subplot(1,2,2), cm.plot('eye movement path','hf',gca);
title(sprintf('Noise %s',cm.noiseFlag));

%% Compute the photocurrent

cm.computeCurrent;

baseFrame   = 25;
baseIMG     = cm.current(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(cm.tSamples,1);

% Compute the error, but because of the ramp up use 25 as the base frame
for ii=1:cm.tSamples
    thisIMG     = cm.current(:,:,ii);
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

%% Get rid of the warm up samples

% Plot the size of the displacement and the error on the same graph
start = baseFrame;
vcNewGraphWin;
plot(d(start:end),eAbsorb(start:end),'o-'); grid on; set(gca,'ylim',[-10 10]);
xlabel('Distance (cones)'); ylabel('Percent error');


%% END