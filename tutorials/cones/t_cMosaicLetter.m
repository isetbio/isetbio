%% Show a letter on a display imaged at the level of the cone mosaic
%
%
% TODO:
%   Write a routine that takes the outline of the letter and superimposes
%   it on the color image of the cm.
%
%   Implement this for the hex mosaic.
%
% Wandell, 2019
%
% See also
%   sceneCreate('letter', ...), displayCreate, ....
%

%%
ieInit;

%% Create a letter on a display

% family, size, dpi
font = fontCreate('A', 'Georgia', 14, 96);
display = 'LCD-Apple';
scene = sceneCreate('letter', font, display);
scene = sceneSet(scene,'wangular',0.6);

% We should pad the scene so the eye movements do not move the scene beyond
% the array

% Here is the scene
% sceneWindow(scene);

%% Push the scene through human optics

oi = oiCreate;
oi = oiCompute(oi,scene);
%oiWindow(oi);
%%  Now image it on the cone mosaic with some fixational eye movements
cm = mosaicLoad('cmosaic_0.5-0.5_0.0-0.0.mat');
% cm = cMosaic('sizeDegs',[0.5, 0.5],'eccentricityDegs',[0,0]);
cm.visualize();

%% Illustrate the cone excitations

cm.integrationTime = 5/1000;  % 5 ms
excitations = cm.compute(oi);

params = cm.visualize('params');
% cm.visualize('help');

params.activation = excitations.^0.5;
params.activationColorMap = hot(1024);
params.verticalActivationColorBar = true;
cm.visualize(params);

%% Add eye movements
eyeMovementDurationSeconds = 100/1000;
cm.emGenSequence(eyeMovementDurationSeconds, ...
        'microsaccadeType', 'none', ...
        'nTrials', 1, ...
        'randomSeed', 10);
    
%% Compute 128 noisy response instances of cone excitation response to the same eye movement path
instancesNum = 128;
% (Instances, Time Samples, Cone index)
[excitationsInstances, noisyExcitationResponseInstances, ~,~,timeAxis] = cm.compute(oi, ...
    'withFixationalEyeMovements', true, ...
    'nTrials', instancesNum);

%% Visualize time-series response of a single cone

% Find the cone with max noise-free response
[~,idx] = max(excitationsInstances(:));
[~,~,targetConeID] = ind2sub(size(excitationsInstances), idx);


ieNewGraphWin;

% Plot the time series response for individual trials
plot(timeAxis, squeeze(noisyExcitationResponseInstances(:,:,targetConeID)), 'k.');
hold on;

% Plot the time series response for the mean of the individual instances
plot(timeAxis, squeeze(mean(noisyExcitationResponseInstances(:,:,targetConeID),1)), 'g-', 'LineWidth', 2.0);
hold on;

% Overlay the noise-free time series response in red
plot(timeAxis, squeeze(excitationsInstances(:,:,targetConeID)), 'r', 'LineWidth', 1.5);
xlabel('time (seconds)');
ylabel('excitations per integration time');
set(gca, 'FontSize', 16);

%% Plot a movie of the excitations

% You can play this video with VLC
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
for ii=1:numel(timeAxis)
    [~,hdl] = cm.plot('excitations',excitationsInstances(1,ii,:));
    drawnow
    thisImg = getframe(gcf);
    writeVideo(vidfile,thisImg);
end
close(vidfile)

%% END

