function varargout = v_osBioPhysObject(varargin)
% Validate the biophysical model of the cone outer segments.
%
% This script tests out our implementation of a biophysical model of the
% isomerizations to photocurrent transduction process that occurs in the
% cone outer segments.
%
% 6/xx/2015    fmr   Created.
% 7/xx/2015    jrg   Test with ISETBIO outersegment object
% 12/31/15     dhb   Added local copy of coneAdapt.            
% 1/7/16       dhb   Rename.  Started to remove reference to coneAdapt.  
%                    Last version with coneAdapt comparison is in tagged
%                    version OSObjectVsOrigValidation.

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;

%% Impulse response in dark
%
% We'll compare the original code and the outer segment object for this
% case.
%
% Set up parameters for stimulus.
nSamples = 2000;        % 2000 samples
timeStep = 1e-4;        % time step
flashIntens = 50000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)

% Create human sensor.
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);

% Create stimulus.
stimulus = zeros(nSamples, 1);
stimulus(1) = flashIntens;
stimulus = reshape(stimulus, [1 1 nSamples]);

% Set photon rates. This is a kluge that appeared
% just for this test, and that should probably go
% away again. This is an artifact of directly specifying the stimulus
% in the sensor, and will not be an issue when the sensor
% is the result of a sensorCompute command on a scene and oi.
sensor = sensorSet(sensor, 'photon rate', stimulus);

% Create outersegment object and get the adapted response.
noiseFlag = 0;
adaptedOS = osBioPhys();
adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag);
adaptedOS = osCompute(adaptedOS, sensor);
osAdaptedCur = osGet(adaptedOS,'coneCurrentSignal');

% Make a plot of what happened using osPlot.
if (runTimeParams.generatePlots)   
    osPlot(adaptedOS, sensor);
end

% Another plot of what happened
if (runTimeParams.generatePlots)
    vcNewGraphWin; hold on
    tme = (1:nSamples)*timeStep;
    plot(tme,squeeze(osAdaptedCur),'r','LineWidth',2);
    xlabel('Time (sec)','FontSize',14);
    ylabel('Photocurrent (pA)','FontSize',14);
    title('Impulse response in the dark','FontSize',16);
end

% Save validation data
tolerance = 0;
UnitTest.validationData('osAdaptedCur',osAdaptedCur);
clear adaptedOS paramsOS

%% Steps + flashes

% Set up parameters for stimulus.
stimPeriod = [4000 22000];          % start and end time of step
nSamples = stimPeriod(2)+4000;      % number of samples
timeStep = 1e-4;                    % time step
flashTime = stimPeriod - 3000;      % time of flashes
flashDur = 10;                      % flash duration (bins)
flashIntens = 10000;                % flash intensity R*/cone/sec

% Create human sensor.
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);

% Set up plot
if (runTimeParams.generatePlots)
    vcNewGraphWin;
    subplot(1, 3, 1); hold on;
    subplot(1, 3, 2); hold on;
    subplot(1, 3, 3); hold on
end

% Go through series of step intensities (slow).
clear SSCur stimIntensity FlashAmp
FlashScFact = 1;                     % scaling of flash to compensate adaptation
nStepIntensities = 11;
for step = 1:nStepIntensities
    
    % Create stimulus: step alone.
    stimIntensity(step) = 50 * 2^step;
    stimulus = zeros(nSamples, 1);
    stimulus(stimPeriod(1):stimPeriod(2)) = stimIntensity(step);
    stimulus = reshape(stimulus, [1 1 nSamples]);
    
    % Increase flash strength if we are in range where cones are adapting.
    if (stimIntensity(step) > 4000);
        FlashScFact = FlashScFact * 2;
    end
    
    % Set photon rates.
    % This is an artifact of directly specifying the stimulus
    % in the sensor, and will not be an issue when the sensor
    % is the result of a sensorCompute command on a scene and oi.
    sensor = sensorSet(sensor, 'photon rate', stimulus);
    
    % Compute using outersegment object osBioPhys.
    noiseFlag = 0;
    adaptedOSStepOnly = osBioPhys();
    adaptedOSStepOnly = osSet(adaptedOSStepOnly, 'noiseFlag', noiseFlag);   
    adaptedOSStepOnly = osCompute(adaptedOSStepOnly, sensor);
    
    % Create stimulus: step + flashes.
    stimulus = zeros(nSamples, 1);
    stimulus(stimPeriod(1):stimPeriod(2)) = stimIntensity(step);
    stimulus(flashTime(1):flashTime(1)+flashDur) = stimulus(flashTime(1):flashTime(1)+flashDur) + flashIntens;
    stimulus(flashTime(2):flashTime(2)+flashDur) = stimulus(flashTime(2):flashTime(2)+flashDur) + flashIntens*FlashScFact;
    stimulus = reshape(stimulus, [1 1 nSamples]);
    
    % Set photon rates.
    sensor = sensorSet(sensor, 'photon rate', stimulus);
    
    % Create outersegment object.
    noiseFlag = 0;
    paramsOS.bgVolts = 0;
    adaptedOS = osBioPhys();
    adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag);
    adaptedOS = adaptedOS.compute(sensor, paramsOS);
    
    % Plot.
    if (runTimeParams.generatePlots)
        subplot(1, 3, 1);
        plot((1:nSamples)*timeStep, adaptedOS.coneCurrentSignal(:), 'r', 'lineWidth', 2);

        subplot(1, 3, 2);
        plot((1:nSamples)*timeStep, adaptedOSStepOnly.coneCurrentSignal(:), 'r', 'lineWidth', 2);

        subplot(1, 3, 3);
        plot((1:nSamples)*timeStep, adaptedOS.coneCurrentSignal(:) - adaptedOSStepOnly.coneCurrentSignal(:), 'r', 'lineWidth', 2);
      
    end
    
    % Summary measures.    
    %tempVar = adaptedCur(:) - adaptedCurStepOnly(:);
    %flashAmp(step) = max(tempVar(flashTime(2):flashTime(2)+1000)) / (FlashScFact * max(tempVar(flashTime(1):flashTime(1)+1000)));
    %ssCur(step) = -(adaptedCurStepOnly(1, 1, stimPeriod(2)) - adaptedCurStepOnly(1, 1, 1))/adaptedCurStepOnly(1, 1, 1);
    %[maxVal, maxLoc] = max(tempVar(flashTime(2):flashTime(2)+1000));
    %tPeak(step) = maxLoc*timeStep;
    
    tempOS = adaptedOS.coneCurrentSignal(:) - adaptedOSStepOnly.coneCurrentSignal(:);
    flashAmpOS(step) = max(tempOS(flashTime(2):flashTime(2)+1000)) / (FlashScFact * max(tempOS(flashTime(1):flashTime(1)+1000)));
    ssCurOS(step) = -(adaptedOSStepOnly.coneCurrentSignal(1, 1, stimPeriod(2)) - adaptedOSStepOnly.coneCurrentSignal(1, 1, 1))/adaptedOSStepOnly.coneCurrentSignal(1, 1, 1);
    
    % Fit transient of step response.
    transientTime = stimPeriod(1)+round(6e-2/timeStep):stimPeriod(2);
    %transient{step} = adaptedCurStepOnly(transientTime);
    %transient{step} = transient{step} - transient{step}(end);
    transientOS{step} = adaptedOSStepOnly.coneCurrentSignal(stimPeriod(1)+round(6e-2/timeStep):stimPeriod(2));
    transientOS{step} = transientOS{step} - transientOS{step}(end);
    
end

% Label plot.
if (runTimeParams.generatePlots)
    subplot(1, 3, 1);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('Adapted Current');
    
    subplot(1, 3, 2);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('Adapted Current Step Only');
    
    subplot(1, 3, 3);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('Adapted Current Steady State');
end

% Compare to model fit.

% Fit current transient - two temporal components which contrast fast and
% slow Calcium feedback terms, 3:1 ratio of amps, rate constants 1 and 12.
%
% This fit not currently being used.
tme = (1:length(transientOS{end}))*timeStep;
coef = [1 1 3 12];
stepcoefOS = nlinfit(tme', squeeze(transientOS{end}), 'sumoftwodblexponentials', coef);
transientFitOS = sumoftwodblexponentials(stepcoefOS,  tme);

% Plot steady-state stimulus-response relation.
% Half max 45000 (Dunn et al 2007)
% Fit steady-state  vs intensity relation for each way of computing
coef = [45000 1];
fitcoefOS = nlinfit(stimIntensity, ssCurOS, 'hill', coef);
fitOS = hill(fitcoefOS, stimIntensity);
if (runTimeParams.generatePlots)
    vcNewGraphWin; 
    hold on;
    semilogx(stimIntensity, ssCurOS, 'ro');
    semilogx(stimIntensity, fitOS,'r');
    xlabel('Background','fontsize',16);
    ylabel('Steady State Current','fontsize',16);
    axis tight;
    title('Steady State Current vs. Background Intensity','fontsize',16);
end

% Plot sensitivity vs background intensity relation
% Half desens around 2500 (Angueyra and Rieke, 2013)
% Fit sensitivity vs intensity relation computed both ways
% coneAdapt data:
coef = [2500];
wfcoefOS = nlinfit(stimIntensity, log10(flashAmpOS), 'weberFechner', coef);
fitOS = weberFechner(wfcoefOS, stimIntensity);
if (runTimeParams.generatePlots)
    vcNewGraphWin; 
    hold on;
    loglog(stimIntensity, flashAmpOS, 'ro');
    loglog(stimIntensity, 10.^fitOS,'r');
    xlabel('Background','fontsize',16);
    ylabel('Sensitivity','fontsize',16);
    axis tight;
    title('Sensitivity vs. Background Intensity','fontsize',16);
end

% Plot transient from last step simulated
if (runTimeParams.generatePlots)
    vcNewGraphWin; hold on
    plot(transientTime, transientOS{end}(:),'r','Linewidth',2);
    xlabel('Time (sec)','fontsize',16);
    ylabel('Photocurrent (pA)','fontsize',16);
    title('Transient for largest step','fontsize',16);
end

% Report fits.  This printout is a little dangerous because the variable
% names are pretty generic and might be viewed as temporary somewhere
% above.
fprintf(1, 'half max step amp = %d\nhalf desensitizing background = %d\nstep time constants = %d and %d\n', fitcoefOS(1), wfcoefOS(1), stepcoefOS(2), stepcoefOS(4));
fprintf(1, 'targets: half max step amp = 45000 R*/cone/sec\n\thalf desensitizing background = 2500 R*/cone/sec \n\tstep time constants = 1 and 12 sec\n');

tolerance = 1e-16;
UnitTest.validationData('osAdaptedCur',osAdaptedCur);
clear adaptedOS adaptedOSSteponly paramsOS paramsOSStepOnly

%% Saccade-like stimuli

% Load experimental data using RDT.
% client = RdtClient(getpref('isetbio','remoteDataToolboxConfig'));
client = RdtClient('isetbio');
client.crp('resources/data/cones');
[eyeMovementExample, eyeMovementExampleArtifact] = client.readArtifact('eyeMovementExample', 'type', 'mat');

% Get mean subtracted current and stimulus from read data
measuredCur = eyeMovementExample.data.Mean;
measuredCur = measuredCur - measuredCur(1,end);
stimulus = eyeMovementExample.data.Stim;
nSamples = length(stimulus);
stimulus = reshape(stimulus, [1 1 nSamples]);

% Compute adapted current using DEQ model.

% Create human sensor.
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', 5e-5);

% Set photon rates.
% This is an artifact of directly specifying the stimulus
% in the sensor, and will not be an issue when the sensor
% is the result of a sensorCompute command on a scene and oi.
sensor = sensorSet(sensor, 'photon rate', stimulus);

% Create outersegment object.
noiseFlag = 0;
adaptedOS = osBioPhys();
adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag);
adaptedOS = osCompute(adaptedOS, sensor);
osAdaptedCur = osGet(adaptedOS, 'coneCurrentSignal');

% Plot the two calculations and compare against measured data.
if (runTimeParams.generatePlots)
    vcNewGraphWin; hold on;
    % Correct for DC offset to measuredCur by setting last value to zero.
    plot((1:nSamples)*5e-5, measuredCur(:)-measuredCur(end), 'b', 'LineWidth', 4);
    plot((1:nSamples)*5e-5, osAdaptedCur(:)-osAdaptedCur(end), 'r', 'LineWidth', 2);
    axis tight
    xlabel('Time (sec)','FontSize',16);
    ylabel('Photocurrent (pA)','FontSize',16);
    title('Comparison for end subtracted eye movement example data','FontSize',14);
    legend('measurements', 'osBioPhys');
end

% Save validation data
tolerance = 0;
UnitTest.validationData('osAdaptedCur',osAdaptedCur);
clear adaptedOS paramsOS

%% Increment/decrement asymmetry

% Set up parameters for stimulus.
stimPeriod = [5000 7500]; % start and end time of step
nSamples = stimPeriod(2)+5000;
timeStep = 2e-4;
Contrast = 1;

% Create human sensor.
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', timeStep);

% Start figure for this next bit
if (runTimeParams.generatePlots)  
    vcNewGraphWin;
end

% Predict responses to increments and decrements across range of light
% levels.
clear MaxInc MaxDec stimIntensity
stepLevels = 7;
for step = 1:stepLevels
    % Create stimulus: step along.
    stimIntensity(step) = 250 * 2^step;          % photons per second
    stimulusInc = zeros(nSamples, 1);
    stimulusInc(100:nSamples-100) = stimIntensity(step);
    stimulusInc(stimPeriod(1):stimPeriod(2)) = stimIntensity(step)*(1 + Contrast);
    stimulus = reshape(stimulusInc, [1 1 nSamples]);
    
    % Set photon rates.
    % This is an artifact of directly specifying the stimulus
    % in the sensor, and will not be an issue when the sensor
    % is the result of a sensorCompute command on a scene and oi.
    sensor = sensorSet(sensor, 'photon rate', stimulus);
    
    % Create outersegment object and compute with it.
    noiseFlag = 0;
    adaptedOSInc = osBioPhys();
    adaptedOSInc = osSet(adaptedOSInc, 'noiseFlag', noiseFlag);    
    adaptedOSInc = osCompute(adaptedOSInc, sensor);    
    osAdaptedCur = osGet(adaptedOSInc, 'coneCurrentSignal');
     
    %  create stimulus: step + flashes
    stimulusDec = zeros(nSamples, 1);
    stimulusDec(100:nSamples-100) = stimIntensity(step);
    stimulusDec(stimPeriod(1):stimPeriod(2)) = stimIntensity(step)*(1 - Contrast);
    stimulus = reshape(stimulusDec, [1 1 nSamples]);
    
    % Set photon rates.
    sensor = sensorSet(sensor, 'photon rate', stimulus);
    
    % And with outersegment object.
    noiseFlag = 0;
    adaptedOSDec = osBioPhys('noiseFlag', noiseFlag);
    paramsOSDec.bgVolts = 0;
    adaptedOSDec = adaptedOSDec.compute(sensor, paramsOSDec);
    
    maxOSInc(step) = adaptedOSInc.coneCurrentSignal(1, 1, stimPeriod(2)-1) - adaptedOSInc.coneCurrentSignal(1, 1, stimPeriod(1)-1);
    maxOSDec(step) = adaptedOSDec.coneCurrentSignal(1, 1, stimPeriod(2)-1) - adaptedOSDec.coneCurrentSignal(1, 1, stimPeriod(1)-1);
    
    % Add current to the figure.
    if (runTimeParams.generatePlots)
        subplot(1,2,2); hold on;
        plot((1:nSamples)*timeStep, adaptedOSInc.coneCurrentSignal(:), (1:nSamples)*timeStep, adaptedOSDec.coneCurrentSignal(:));
    end
end

% Finish off the figure
if (runTimeParams.generatePlots)
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('osBioPhys');
end

if (runTimeParams.generatePlots)
    vcNewGraphWin; 
    hold on;
    semilogx(stimIntensity, -maxOSDec ./ maxOSInc, 'ro');
    xlabel('Background','fontsize',16);
    ylabel('Dec/Inc Response Ratio','fontsize',16);
    title('Dec/Inc Response Ratio vs. Background Intensity','fontsize',14);
end

% Tuck data away.
tolerance = 1e-12;
UnitTest.validationData('osAdaptedCurInc',adaptedOSInc.coneCurrentSignal);
UnitTest.validationData('osAdaptedCurDec',adaptedOSDec.coneCurrentSignal);

end

