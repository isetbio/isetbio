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

hparam = harmonicP('freq',[1, 7],'ang',[0 0],...
    'ph',[pi/2, pi/2], 'contrast',[1 1],...
    'row',256,'col',256);
scene = sceneCreate('harmonic',hparam);
scene = sceneSet(scene,'fov',0.6);
sceneWindow(scene);

%% Push the scene through human optics

oi = oiCreate('human wvf');
oi = oiCompute(oi,scene,'pad value','mean');
%oiWindow(oi);

%%  Now image it on the cone mosaic with some fixational eye movements

% cm = mosaicLoad('cmosaic_0.5-0.5_0.0-0.0.mat'); ecc = 0;
% cm = mosaicLoad('cmosaic_0.5-0.5_10.0-0.0.mat'); ecc = 10;

cm = mosaicLoad('cmosaic_1.0-1.0_0.0-0.0.mat'); ecc = 0;
% cm = mosaicLoad('cmosaic_1.0-1.0_10.0-0.0.mat'); ecc = 10;
cm.visualize();

%% Illustrate the cone excitations

cm.integrationTime = 5/1000;  % 5 ms
[~,excitations] = cm.compute(oi);

params = cm.visualize('params');
% cm.visualize('help');

params.activation = excitations.^0.5;
params.activationColorMap = hot(1024);
params.verticalActivationColorBar = true;
cm.visualize(params);

%% Show a plot
cm.plot('excitations',excitations);
cm.plot('excitations horizontal line',excitations,'ydeg',0,'cone type','L');
%% Add eye movements
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

%% Plot a movie of the excitations
fname = sprintf('%s-harmonic-ecc-%d.mp4',fullfile(isetRootPath,'local','excitations'),ecc);
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

fname = sprintf('%s-harmonic-ecc-%d.mp4',fullfile(isetRootPath,'local','current'),ecc);
cm.movie(timeAxis,current,'filename',fname);
implay(fname);

%% END

