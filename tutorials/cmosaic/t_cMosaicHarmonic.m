%% Show harmonics on the cone mosaic
%
% Wandell, 2019
%
% See also
%

%%
ieInit;

%% Create the sume of two harmonics on a display
% The first harmonic is a low frequency horizontal grating.  The second is a high frequency horizontal grating.  The two are summed to create a visual impression of high contrast and low contrast from left to right.

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

cm = mosaicLoad('cmosaic_0.5-0.5_0.0-0.0.mat'); ecc = 0;
%
%   cm = mosaicLoad('cmosaic_0.5-0.5_10.0-0.0.mat'); ecc = 10;
%
cm.plot('mosaic');

%% Illustrate the cone excitations

cm.integrationTime = 5/1000;  % 5 ms
[~,excitations] = cm.compute(oi);

cm.plot('excitations', excitations, ...
    'plotTitle', 'Cone excitations (static)');

%% Show a plot of the horizontal line profile
cm.plot('excitations horizontal line',excitations,...
    'ydeg',0,...
    'cone type','L');

%% Add eye movements
% We use a short duration (100 ms) so the tutorial runs quickly.
eyeMovementDurationSeconds = 100/1000;
cm.emGenSequence(eyeMovementDurationSeconds, ...
    'microsaccadeType', 'none', ...
    'nTrials', 1, ...
    'randomSeed', 10);

%% Compute noisy cone excitations to the same eye movement path
instancesNum = 1;
% (Instances, Time Samples, Cone index)
[~, excitations, ~,~,timeAxis] = cm.compute(oi, ...
    'withFixationalEyeMovements', true, ...
    'nTrials', instancesNum);

%% Visualize time-series response of a single cone
% Find the cone with the maximum mean excitation across time and trials
[~, targetConeID] = max(mean(excitations, [1 2]));

% Plot the time series for that cone. Pass coneindices to select it.
cm.plot('time series', excitations, ...
    'timeaxis', timeAxis, ...
    'coneindices', targetConeID);

%% Plot a movie of the excitations
% To ensure the movie writes properly and runs quickly, we use a smaller
% mosaic and fewer frames in this tutorial.
fname = sprintf('%s-harmonic-ecc-%d.mp4',fullfile(isetRootPath,'local','excitations'),ecc);

% Generating the movie can take a moment for large mosaics
cm.movie(timeAxis,excitations,'filename',fname);

% Play the movie (unless we are publishing the script)
stackInfo = dbstack;
if ~any(strcmp({stackInfo.name}, 'publish'))
    implay(fname);
end

%%  Photocurrent -> Compute the photocurrent from the excitations

meanIso = cm.meanIsomerizations(excitations);

% Mean isomerizations per temporal sampling bin
irf = currentIRF(meanIso*cm.integrationTime,ecc,timeAxis);

% Convolution.
current = cm.current(excitations,irf,timeAxis);

fname = sprintf('%s-harmonic-ecc-%d.mp4',fullfile(isetRootPath,'local','current'),ecc);
cm.movie(timeAxis,current,'filename',fname);

% Play the movie (unless we are publishing the script)
stackInfo = dbstack;
if ~any(strcmp({stackInfo.name}, 'publish'))
    implay(fname);
end


