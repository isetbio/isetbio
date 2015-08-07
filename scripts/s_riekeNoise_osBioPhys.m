%% s_coneAdaptNoise
%  This script computes and analyzes noise in cone adaptation. 
%     
%  The stimulus used in the experiment is a step function (dark - light
%     - dark). The measurement data is in NoiseData.mat
%
%  The cone adapted current is computed with the physiological differential
%  equation model (see riekeAdapt) and noise is computed as the difference
%  of the measurements and estimates.
%
%  Noise autocorrelation function and spectral density are computed
%
%  Note:
%    The coneAdapte(sensor, 'rieke') cannot predict the response currently
%    HJ will check and get some new numbers from Fred.
%
%  (HJ) ISETBIO TEAM, 2014

%% Init
s_initISET;

%% Generate stimulus
%  setup parameters for stimulus
nSamples      = 140000;         % total samples in the stimulus
stimPeriod    = [40001 120000]; % start and end time of the light
stimIntensity = 78435;          % photons per second

%  create stimulus
stimulus = zeros(nSamples, 1);
stimulus(stimPeriod(1):stimPeriod(2)) = stimIntensity;
stimulus = reshape(stimulus, [1 1 nSamples]);

% plot stimulus
vcNewGraphWin; plot(stimulus(:));

%% Load measurement data
pathNoiseData = '/home/james/Downloads/fred cone model/';
noiseData = load([pathNoiseData 'NoiseData.mat']);
% load('NoiseData.mat');
measuredCur = noiseData.StepsExample.Data{3}; % current in pA

%% Compute adapted current with differential equations
%  create human sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone

% set photon rates
sensor = sensorSet(sensor, 'photon rate', stimulus);

% compute adapted current
params.bgVolts  = 0;
[~, adaptedCur] = coneAdapt(sensor, 'rieke', params);
% adaptedCur = adaptedCur + (measuredCur(1)-adaptedCur(1));
% plot against measurement data
vcNewGraphWin; hold on;
plot(1:nSamples, adaptedCur(:), '--r', 'lineWidth', 2);
plot(1:nSamples, measuredCur(1,:));

vcNewGraphWin; plot(adaptedCur(:), measuredCur(1,:), '.');

% Create outersegment object
noiseFlag = 0;
adaptedOS = osBioPhys('noiseFlag', noiseFlag);
% adaptedOS = osLinear('noiseflag', noiseFlag);
paramsOS.bgVolts = params.bgVolts;
% paramsOS.dc = 0; % removes dc
adaptedOS = adaptedOS.compute(sensor, paramsOS);
% adaptedOS.plotResults(sensor)

vcNewGraphWin; plot(squeeze(adaptedOS.ConeCurrentSignal), measuredCur(1,:), '.');

% figure; scatter(adaptedOS.ConeCurrentSignal,adaptedCur)

%% Comptue and analyze noise
%  Compute noise
noise = measuredCur(1,:) - squeeze(adaptedCur)';

%  Plot ACF and PACF
vcNewGraphWin; 
subplot(1, 2, 1); autocorr(noise);
subplot(1, 2, 2); parcorr(noise);

% Plot power spectral density

% Fit ARIMA model to the noise
mdl = arima(3, 0, 2);
mdl = estimate(mdl, noise);

res = infer(mdl, noise);

% see how good the fit is
vcNewGraphWin;
subplot(2,2,1); plot(res ./ sqrt(mdl.Variance)); 
title('standardized residual');
subplot(2,2,2);  qqplot(res);   % check normality, should fall along a line
subplot(2,2,3);  autocorr(res); % independence, should be small
subplot(2,2,4);  parcorr(res);  % independence, should be small

%  Compute noise for OS object
noiseOS = measuredCur(1,:) - squeeze(adaptedOS.ConeCurrentSignal)';

%  Plot ACF and PACF
vcNewGraphWin; 
subplot(1, 2, 1); autocorr(noise);
subplot(1, 2, 2); parcorr(noise);

% Plot power spectral density

% Fit ARIMA model to the noise
mdlOS = arima(3, 0, 2);
mdlOS = estimate(mdlOS, noiseOS);

resOS = infer(mdlOS, noiseOS);

% see how good the fit is
vcNewGraphWin;
subplot(2,2,1); plot(resOS ./ sqrt(mdlOS.Variance)); 
title('standardized residual');
subplot(2,2,2);  qqplot(resOS);   % check normality, should fall along a line
subplot(2,2,3);  autocorr(resOS); % independence, should be small
subplot(2,2,4);  parcorr(resOS);  % independence, should be small

%% Generate noise to simulation data
%  add noise
paramsNoise.sampTime = 1e5;
adaptedCurNoise = riekeAddNoise(adaptedCur,params)*1;

% Create outersegment object
clear adaptedOS
noiseFlag = 1;
adaptedOS = osBioPhys('noiseFlag', noiseFlag);
% adaptedOS = osLinear('noiseflag', noiseFlag);
paramsOS.bgVolts = params.bgVolts;
% paramsOS.dc = 0; % removes dc
adaptedOS = adaptedOS.compute(sensor, paramsOS);

vcNewGraphWin; hold on;
plot(1:nSamples, adaptedCur(:), '--r', 'lineWidth', 2);
plot(1:nSamples, squeeze(adaptedOS.ConeCurrentSignalPlusNoise), '--k', 'lineWidth', 2);
plot(1:nSamples, measuredCur(1,:));

plot(1:nSamples, adaptedCurNoise(:), '--g', 'lineWidth', 2);

% figure;
% scatter(adaptedCur(:),squeeze(adaptedOS.ConeCurrentSignalPlusNoise));
%% END