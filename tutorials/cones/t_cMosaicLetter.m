%% Simulate a letter imaged at the level of the cone mosaic
%
% The letter is simulated on a display. It is then rendered through the
% optics onto the modern cone mosaic (cMosaic) object.
%
% See also
%   sceneCreate('letter', ...), displayCreate, ....
%

% Things TODO:
%   Write a routine that takes the outline of the letter and superimposes
%   it on the color image of the cm.
%
%   Implement this for the hex mosaic.
%

%%
ieInit;

%% Create a letter on a display

% family, size, dpi
font = fontCreate('A', 'Georgia', 18, 96);
display = 'LCD-Apple';
scene = sceneCreate('letter', font, display);
scene = sceneSet(scene,'wangular',0.5);

%{
% To make a second letter and combine, perhaps?

% family, size, dpi
scene = sceneCrop(scene,[96 1 sz(2)-97 sz(1)-1]);
sz = sceneGet(scene,'size');

font = fontCreate('B', 'Georgia', 18, 96);
display = 'LCD-Apple';
scene2 = sceneCreate('letter', font, display);
scene2 = sceneSet(scene2,'wangular',0.5);
scene2 = sceneCrop(scene2,[96 1 sz(2)-65 sz(1)-1]);
scene = sceneCombine(scene,scene2,'direction','horizontal');
%}

sceneWindow(scene);

%% Calculate through human optics

oi = oiCreate;
oi = oiCompute(oi,scene);

%oiWindow(oi);

%%  Now image it on the cone mosaic with some fixational eye movements

% Here is a stored cMosaic
cm = mosaicLoad('cmosaic_0.5-0.5_0.0-0.0.mat'); ecc = 0;
% cm = mosaicLoad('cmosaic_0.5-0.5_10.0-0.0.mat'); ecc = 10;
% cm = mosaicLoad('cmosaic_1.0-1.0_0.0-0.0.mat'); ecc = 0;
% cm = mosaicLoad('cmosaic_1.0-1.0_10.0-0.0.mat'); ecc = 10;

% Visualize the loaded mosaic
cm.visualize();

%% Calculate the cone excitations from a brief exposure

cm.integrationTime = 5/1000;  % 5 ms
[~,excitations] = cm.compute(oi);

% To control the visualization, we set these parameters
% cm.visualize('help');

params = cm.visualize('params');
params.activation = excitations.^0.5;
params.activationColorMap = hot(1024);
params.verticalActivationColorBar = true;

cm.visualize(params);

%% Visualize the impact of eye movements

% These calculations can take a little while (few minutes)
eyeMovementDurationSeconds = 400/1000;
cm.emGenSequence(eyeMovementDurationSeconds, ...
        'microsaccadeType', 'none', ...
        'nTrials', 1, ...
        'randomSeed', 10);
    
%% Compute noisy response instances of cone excitation response to the same eye movement path

instancesNum = 1;

% (Instances, Time Samples, Cone index)
[~, excitations, ~,~,timeAxis] = cm.compute(oi, ...
    'withFixationalEyeMovements', true, ...
    'nTrials', instancesNum);

%% Plot a movie of the excitations
fname = sprintf('%s-ecc-%d.mp4',fullfile(isetRootPath,'local','excitations'),ecc);
cm.movie(timeAxis,excitations,'filename',fname);
implay(fname);

% TODO:  Figure this out.
%
% cmd = sprintf('vlc %s',mp4File);
% [s,r] = system(cmd)
%
% implay(mp4File);
% You can also import the file into PowerPoint or use VLC
%

%%  Photocurrent

% For the coneMosaic we used to be able to do this.
%
%  thisOS = osBioPhys;
%  current = thisOS.osCompute(cMosaic);
%
meanIso = cm.meanIsomerizations(excitations);

% Mean isomerizations per temporal sampling bin
irf = currentIRF(meanIso*cm.integrationTime,ecc,timeAxis);

% Convolution.
current = cm.current(excitations,irf,timeAxis);

fname = sprintf('%s-ecc-%d.mp4',fullfile(isetRootPath,'local','current'),ecc);
cm.movie(timeAxis,current,'filename',fname);
implay(fname);

%% END

%% Visualize time-series response of a single cone
%{
% Find the cone with max noise-free response
[~,idx] = max(excitations(:));
[~,~,targetConeID] = ind2sub(size(excitations), idx);

ieNewGraphWin;

% Plot the time series response for individual trials
plot(timeAxis, squeeze(excitations(:,:,targetConeID)), 'b--');
hold on;

% Plot the time series response for the mean of the individual instances
plot(timeAxis, squeeze(mean(excitations(:,:,targetConeID),1)), 'k-', 'LineWidth', 2.0);
hold on;

xlabel('Time (seconds)');
ylabel('Excitations per integration time');
set(gca, 'FontSize', 16);
grid on;
%}

