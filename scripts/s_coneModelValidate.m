%% s_coneModelValidate
% 
% This script compares the DEQ based cone model output with several key
% measured properties of the cone responses:
%   - impulse response in the dark
%   - responses to steps + flashes to measure kinetics of step response and
%   adaptation of step and flash responses
%   - responses to stimuli designed to mimic those encountered by a cone
%   during saccades about a natural scene
%   - responses to light increments and decrements
%
%   references for original data:
%       Dunn et al. (2007)
%       Angueyra and Rieke (2013)
%
% Created 6/2015 FMR



%% Init
s_initISET;

%% impulse response in dark
%  set up parameters for stimulus
nSamples = 2000;        % 2000 samples
timeStep = 1e-4;        % time step
flashIntens = 50000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)

%  create human sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);

%  create stimulus
stimulus = zeros(nSamples, 1);
stimulus(1) = flashIntens;
stimulus = reshape(stimulus, [1 1 nSamples]);

% set photon rates
sensor = sensorSet(sensor, 'photon rate', stimulus);

% compute model current and baseline correct
params.bgVolts  = 0;
[~, adaptedCur] = coneAdapt(sensor, 'rieke', params);
adaptedCur(:) = adaptedCur(:) - adaptedCur(:, :, nSamples);


% compute fit to measured response (Angueyra and Rieke, 2013)
coef = [1 0.05 0.1 1 0];
tme = (1:nSamples)*timeStep;
impcoef = nlinfit(tme', squeeze(adaptedCur), 'ConeEmpiricalDimFlash', coef);
fit = ConeEmpiricalDimFlash(impcoef, tme');
expcoef = [5 0.02 0.03 0.53 34];            % fit to measured response
expfit = ConeEmpiricalDimFlash(expcoef, tme');

% compare model vs empirical fit to data
figure(1); clf;
tme = (1:nSamples)*timeStep;
plot(tme, fit, tme, squeeze(adaptedCur), tme, expfit);
xlabel('sec');
ylabel('pA');
title('impulse response');
legend('model fit', 'model resp', 'exp fit');

%% steps + flashes

%  set up parameters for stimulus
stimPeriod = [4000 22000];          % start and end time of step
nSamples = stimPeriod(2)+4000;      % number of samples
timeStep = 1e-4;                    % time step
flashTime = stimPeriod - 3000;      % time of flashes
flashDur = 10;                      % flash duration (bins)
flashIntens = 10000;                % flash intensity R*/cone/sec
%  create human sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);

figure(1); clf; subplot(1, 3, 1); hold on; subplot(1, 3, 2); hold on; subplot(1, 3, 3); hold on
xlabel('sec');
ylabel('pA');

clear SSCur stimIntensity FlashAmp
FlashScFact = 1;                     % scaling of flash to compensate adaptation

% go through series of step intensities (slow)
for step = 1:11
    
    %  create stimulus: step alone
    stimIntensity(step) = 50 * 2^step;          
    stimulus = zeros(nSamples, 1);
    stimulus(stimPeriod(1):stimPeriod(2)) = stimIntensity(step);
    stimulus = reshape(stimulus, [1 1 nSamples]);
    
    % increase flash strength if we are in range where cones adapting
    if (stimIntensity(step) > 4000);
        FlashScFact = FlashScFact * 2;
    end
    
    % set photon rates
    sensor = sensorSet(sensor, 'photon rate', stimulus);

    % compute adapted current
    params.bgVolts  = 0;
    [~, adaptedCurStepOnly] = coneAdapt(sensor, 'rieke', params);

    %  create stimulus: step + flashes
    stimulus = zeros(nSamples, 1);
    stimulus(stimPeriod(1):stimPeriod(2)) = stimIntensity(step);
    stimulus(flashTime(1):flashTime(1)+flashDur) = stimulus(flashTime(1):flashTime(1)+flashDur) + flashIntens;
    stimulus(flashTime(2):flashTime(2)+flashDur) = stimulus(flashTime(2):flashTime(2)+flashDur) + flashIntens*FlashScFact;
    stimulus = reshape(stimulus, [1 1 nSamples]);

    % set photon rates
    sensor = sensorSet(sensor, 'photon rate', stimulus);

    % compute adapted current
    params.bgVolts  = 0;
    [~, adaptedCur] = coneAdapt(sensor, 'rieke', params);

    % plot
    subplot(1, 3, 1);
    plot((1:nSamples)*timeStep, adaptedCur(:), 'k', 'lineWidth', 2);
    subplot(1, 3, 2);
    plot((1:nSamples)*timeStep, adaptedCurStepOnly(:), 'k', 'lineWidth', 2);
    subplot(1, 3, 3);
    plot((1:nSamples)*timeStep, adaptedCur(:) - adaptedCurStepOnly(:), 'k', 'lineWidth', 2);
    temp = adaptedCur(:) - adaptedCurStepOnly(:);
    pause(1);
    
    % summary statistics
    FlashAmp(step) = max(temp(flashTime(2):flashTime(2)+1000)) / (FlashScFact * max(temp(flashTime(1):flashTime(1)+1000)));
    SSCur(step) = -(adaptedCurStepOnly(1, 1, stimPeriod(2)) - adaptedCurStepOnly(1, 1, 1))/adaptedCurStepOnly(1, 1, 1);
    [MaxVal, MaxLoc] = max(temp(flashTime(2):flashTime(2)+1000));
    TPeak(step) = MaxLoc*timeStep;
    
    % fit transient of step response
    Transient = adaptedCurStepOnly(stimPeriod(1)+round(6e-2/timeStep):stimPeriod(2));
    Transient = Transient - Transient(length(Transient));
    
end

% fit steady-state stimulus-response relation
% half max 45000 (Dunn et al 2007)
figure(2);clf;
semilogx(stimIntensity, SSCur, 'o');
hold on;
coef = [45000 1];
fitcoef = nlinfit(stimIntensity, SSCur, 'hill', coef);
fit = hill(fitcoef, stimIntensity);
semilogx(stimIntensity, fit);
xlabel('background');
ylabel('steady-state current');
axis tight;

% fit sensitivity vs intensity relation
% half desens around 2500 (Angueyra and Rieke, 2013)
figure(3);clf;
loglog(stimIntensity, FlashAmp, 'o');
hold on;
coef = [2500];
wfcoef = nlinfit(stimIntensity, log10(FlashAmp), 'weber_fechner', coef);
fit = weber_fechner(wfcoef, stimIntensity);
loglog(stimIntensity, 10.^fit);
xlabel('background');
ylabel('sensitivity');
axis tight;

% fit current transient - two temporal components which contrast fast and
% slow Calcium feedback terms
% 3:1 ratio of amps, rate constants 1 and 12
tme = (1:length(Transient))*timeStep;
coef = [1 1 3 12];
stepcoef = nlinfit(tme', squeeze(Transient), 'dblexponentialnomean', coef);
fit = dblexponentialnomean(stepcoef, tme);
figure(4); clf
plot(tme, Transient(:), tme, fit);
xlabel('time');
ylabel('pA');

fprintf(1, 'half max step amp = %d\nhalf desensitizing background = %d\nstep time constants = %d and %d\n', fitcoef(1), wfcoef(1), stepcoef(2), stepcoef(4));
fprintf(1, 'targets: half max step amp = 45000 R*/cone/sec\n\thalf desensitizing background = 2500 R*/cone/sec \n\tstep time constants = 1 and 12 sec\n');
%% saccade-like stimuli

% load experimental data and baseline correct
load('~/matlab/isetbio/scripts/ConeData/EyeMovementsExample');
measuredCur = EyeMovementsExample.Mean - EyeMovementsExample.Mean(length(EyeMovementsExample.Mean));

%  create stimulus
stimulus = EyeMovementsExample.Stim;
nSamples = length(stimulus);
stimulus = reshape(stimulus, [1 1 nSamples]);

% Compute adapted current using DEQ model

%  create human sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', 5e-5);

% set photon rates
sensor = sensorSet(sensor, 'photon rate', stimulus);

% compute adapted current
params.bgVolts  = 0;
[~, adaptedCur] = coneAdapt(sensor, 'rieke', params);
adaptedCur = adaptedCur - adaptedCur(:, :, length(stimulus));

% plot against measurement data
figure(5); clf; hold on;
plot((1:nSamples)*5e-5, measuredCur(1,:));
plot((1:nSamples)*5e-5, adaptedCur(:), 'k', 'LineWidth', 2);
axis tight
xlabel('sec');
ylabel('pA');

%% increment/decrement asymmetry

%  set up parameters for stimulus
stimPeriod = [5000 7500]; % start and end time of step
nSamples = stimPeriod(2)+5000;
timeStep = 2e-4;
Contrast = 1;

%  create human sensor
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);

figure(1); clf;
clear MaxInc MaxDec stimIntensity;

% predict responses to increments and decrements across range of light
% levels
for step = 1:7
    
    %  create stimulus: step alone
    stimIntensity(step) = 250 * 2^step;          % photons per second
    stimulusInc = zeros(nSamples, 1);
    stimulusInc(100:nSamples-100) = stimIntensity(step);
    stimulusInc(stimPeriod(1):stimPeriod(2)) = stimIntensity(step)*(1 + Contrast);
    stimulus = reshape(stimulusInc, [1 1 nSamples]);

    % set photon rates
    sensor = sensorSet(sensor, 'photon rate', stimulus);

    % compute adapted current
    params.bgVolts  = 0;
    [~, adaptedCurInc] = coneAdapt(sensor, 'rieke', params);

    %  create stimulus: step + flashes
    stimulusDec = zeros(nSamples, 1);
    stimulusDec(100:nSamples-100) = stimIntensity(step);
    stimulusDec(stimPeriod(1):stimPeriod(2)) = stimIntensity(step)*(1 - Contrast);
    stimulus = reshape(stimulusDec, [1 1 nSamples]);

    % set photon rates
    sensor = sensorSet(sensor, 'photon rate', stimulus);

    % compute adapted current
    params.bgVolts  = 0;
    [~, adaptedCurDec] = coneAdapt(sensor, 'rieke', params);

    % summary stats
    MaxInc(step) = adaptedCurInc(1, 1, stimPeriod(2)-1) - adaptedCurInc(1, 1, stimPeriod(1)-1);
    MaxDec(step) = adaptedCurDec(1, 1, stimPeriod(2)-1) - adaptedCurDec(1, 1, stimPeriod(1)-1);

    % plot
    plot((1:nSamples)*timeStep, adaptedCurInc(:), (1:nSamples)*timeStep, adaptedCurDec(:));
    pause(0.1);
end

figure(2); clf
semilogx(stimIntensity, -MaxDec ./ MaxInc, 'o');
xlabel('background');
ylabel('dec/inc');
