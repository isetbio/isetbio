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
[~,excitations] = cm.compute(oi);

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
[~, excitations, current,~,timeAxis] = cm.compute(oi, ...
    'withFixationalEyeMovements', true, ...
    'nTrials', instancesNum);

%% Visualize time-series response of a single cone

% Find the cone with max noise-free response
[~,idx] = max(excitations(:));
[~,~,targetConeID] = ind2sub(size(excitations), idx);


ieNewGraphWin;

% Plot the time series response for individual trials
plot(timeAxis, squeeze(excitations(:,:,targetConeID)), 'k--');
hold on;

% Plot the time series response for the mean of the individual instances
plot(timeAxis, squeeze(mean(excitations(:,:,targetConeID),1)), 'g-', 'LineWidth', 2.0);
hold on;

% Overlay the noise-free time series response in red
plot(timeAxis, squeeze(excitationsInstances(:,:,targetConeID)), 'r', 'LineWidth', 1.5);
xlabel('time (seconds)');
ylabel('excitations per integration time');
set(gca, 'FontSize', 16);

%% Plot a movie of the excitations

mp4File = cm.movie(timeAxis,excitations);

%
% implay(mp4File);
% You can also import the file into PowerPoint or use VLC
%
% For the coneMosaic we used to be able to do this.
%
%  thisOS = osBioPhys;
%  current = thisOS.osCompute(cMosaic);
%
% Not sure what NC intends for the cMosaic class.

%% END

% {
  foo = squeeze(mean(excitations(:,:,cm.lConeIndices),1));
  L = mean(foo(:))/cm.integrationTime 
  foo = squeeze(mean(excitations(:,:,cm.mConeIndices),1));
  M = mean(foo(:))/cm.integrationTime 
  foo = squeeze(mean(excitations(:,:,cm.sConeIndices),1));
  S = mean(foo(:))/cm.integrationTime   

% 

cMosaic = coneMosaic('os', osLinear, 'pattern', [2 2 2]);
cMosaic.integrationTime = cm.integrationTime;

%% Setup background levels
% Mean of isomerization stimulus in R*/sec (scaled by time step to be
% placed in bins of photons). The mean rate affects the magnitude of the
% impulse response due to adaptation.
%
% For this tutorial, background rates are the same for all classes of cone.
% Results for the L cones are plotted.
meanIsoArray = [L M S] * cMosaic.integrationTime;

%% Allocate space for the impulse responses
% A litte ugly, we get an impulse response and find its size.
nTimeSamples = size(cMosaic.os.linearFilters(cMosaic), 1);
fovea = zeros(nTimeSamples, length(meanIsoArray));
periphery = zeros(nTimeSamples, length(meanIsoArray));

%%  Loop on different background rates and and compute
os = cMosaic.os;
for ii = 1:length(meanIsoArray)
    % Set mean background for start of current output
    cMosaic.absorptions = repmat(meanIsoArray(ii), 1, 3);

    % Compute outer segment currents for the fovea. Do this by pulling out
    % the first column of return, corresponding to the L cones.
    tmp = os.linearFilters(cMosaic, 'eccentricity', 0);
    fovea(:, ii) = tmp(:, 1);

    % Do it for the periphery
    tmp = os.linearFilters(cMosaic, 'eccentricity', 15);
    periphery(:, ii) = tmp(:, 1);

end

timeAxis = cMosaic.os.timeAxis;

ieNewGraphWin;
plot(timeAxis,fovea);
xlabel('Time (s)');
ylabel('A.U.');
grid on;

ieNewGraphWin;
plot(timeAxis,periphery);
xlabel('Time (s)');
ylabel('A.U.');
grid on;

%}