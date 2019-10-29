%% ----------- Dead leaves calculation ----------------
%
% Experiments with photocurrent and with how we calculate the change from
% frame-to-frame.  Here we use the prior frame as the base, rather than a
% fixed frame at the beginning of the whole period.
%
% We might compare these numbers with what happens during a saccade or
% during a change in the stimulus.  
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
cm.computeCurrent;

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

%% Demonstrate the quadrature shifting (again)

baseFrame   = 2;
baseIMG     = cm.absorptions(:,:,baseFrame);
eBase = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
eAbsorb = zeros(cm.tSamples,1);

% Try just shifting the base frame.
% We get a very solid result.  No noise and pure shifting has no impact
% This was so repeatable, I commented it out.
for ii=1:cm.tSamples
    thisIMG = circshift(baseIMG,ii,2);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% {
% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
plot(d,eAbsorb,'o-'); grid on; set(gca, 'ylim',[-1 1]); 
title('Circular shift')
xlabel('Distance (cones)'); ylabel('Percent error');
drawnow
%}

%% Now use the absorptions and eye movements

% In these calculations, we are comparing the prior 5 msec with the current
% 5 msec of absorptions

for ii=2:cm.tSamples
    baseIMG     = cm.absorptions(:,:,ii-1);
    eBase       = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
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
title('Quadratic case, absorptions, 1-back')

subplot(1,2,2), cm.plot('eye movement path','hf',gca);
title(sprintf('Noise %s',cm.noiseFlag));
drawnow

%% Use the photocurrent response rather than absorptions

startFrame   = 25;   % Eliminate the warm up period
eAbsorb = zeros(cm.tSamples,1);

% Compute the error, but because of the ramp up use 25 as the base frame
for ii=startFrame:cm.tSamples
    baseIMG     = cm.current(:,:,ii-1);
    eBase       = dot(baseIMG(:),sQuad(:))^2 + dot(baseIMG(:),cQuad(:))^2;
    
    thisIMG     = cm.current(:,:,ii);
    eAbsorb(ii) = dot(thisIMG(:),sQuad(:))^2 + dot(thisIMG(:),cQuad(:))^2;
    eAbsorb(ii) = 100*(eAbsorb(ii) - eBase)/eBase;
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
subplot(1,2,1), plot(d(startFrame:end),eAbsorb(startFrame:end),'o-'); 
grid on; set(gca,'ylim',[-10 10]);
xlabel('Distance (cones)'); ylabel('Percent error');
title('Quadratic case, current, 1-back')

subplot(1,2,2), cm.plot('eye movement path','hf',gca);
title(sprintf('Noise %s',cm.noiseFlag));
drawnow

%% Suppose we used the linear responses from the photocurrent

startFrame   = 25;   % Eliminate the warm up period
eAbsorbS = zeros(cm.tSamples,1);
eAbsorbC = zeros(cm.tSamples,1);

% Compute the error, but because of the ramp up use 25 as the base frame
for ii=startFrame:cm.tSamples
    baseIMG     = cm.current(:,:,ii-1);
    thisIMG     = cm.current(:,:,ii);

    eBase        = dot(baseIMG(:),cQuad(:));    
    eAbsorbC(ii) = dot(thisIMG(:),cQuad(:));
    eAbsorbC(ii) = 100*(eAbsorbC(ii) - eBase)/eBase;
    
    eBase        = dot(baseIMG(:),sQuad(:));
    eAbsorbS(ii) = dot(thisIMG(:),sQuad(:));
    eAbsorbS(ii) = 100*(eAbsorbS(ii) - eBase)/eBase;
    
    % fprintf('Difference (percentage) %f\n',eAbsorb(ii));
end

% Plot the size of the displacement and the error on the same graph
d = sqrt(cm.emPositions(:,1).^2 + cm.emPositions(:,2).^2);
vcNewGraphWin;
subplot(1,2,1), 
plot(d(startFrame:end),eAbsorbS(startFrame:end),'ob-', ...
    d(startFrame:end),eAbsorbC(startFrame:end),'xr-'); 
grid on; 
xlabel('Distance (cones)'); ylabel('Percent error');
title('Linear case')
subplot(1,2,2), cm.plot('eye movement path','hf',gca);
title(sprintf('Noise %s',cm.noiseFlag));

%% END