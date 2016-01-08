function varargout = v_osLinearObject(varargin)
%
% Validate the linear model of the cone outer segments.
%
% This script tests out our implementation of a linear model of the
% isomerizations to photocurrent transduction process that occurs in the
% cone outer segments.
%
% 6/xx/2015    fmr   Created.
% 7/xx/2015    jrg   Test with ISETBIO outersegment object
% 12/31/15     dhb   Cleaned up.
%              dhb   Added local copy of coneAdapt and riekeLinearCone.
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
% We'll compare the original code and the outer segment object for this case.
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

% Set photon rates.  This is a kluge that appeared
% just for this test, and that should probably go
% away again. 
% This is an artifact of directly specifying the stimulus
% in the sensor, and will not be an issue when the sensor
% is the result of a sensorCompute command on a scene and oi.
sensor = sensorSet(sensor, 'photon rate', stimulus);

% Compute model current, using the old code.
params.bgVolts  = 0;
[~, adaptedCur] = coneAdapt(sensor, 'linear', params);

% Create outersegment object and get the adapted response.
noiseFlag = 0; 
adaptedOS = osLinear();
adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag);

% Compute the linear response and get the current.
adaptedOS = osCompute(adaptedOS, sensor);
osAdaptedCur = osGet(adaptedOS, 'coneCurrentSignal');

% Make a plot of what happened using osPlot.
if (runTimeParams.generatePlots)   
    osPlot(adaptedOS, sensor);
end

% Plot a comparison of the two. These should be identical since we think
% they are the same code doing the same thing.
if (runTimeParams.generatePlots)
    vcNewGraphWin; hold on
    tme = (1:nSamples)*timeStep;
    plot(tme,squeeze(adaptedCur),'r','LineWidth',4);
    plot(tme,squeeze(osAdaptedCur),'k:','LineWidth',2);
    xlabel('Time (sec)','FontSize',14);
    ylabel('Current (pA)','FontSize',14);
    title('impulse response in the dark','FontSize',16);
    legend('original code model resp', 'osLinear object resp');
end
tolerance = 1e-12;
UnitTest.assertIsZero(max(abs(adaptedCur-osAdaptedCur)),'Comparison for dark impulse response',tolerance);
UnitTest.validationData('adaptedCur',adaptedCur);
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
sensor = sensorSet(sensor,'cone type', 2); % set s cone

% Set up plot
if (runTimeParams.generatePlots)   
    vcNewGraphWin;
    subplot(1, 3, 1); hold on;
    subplot(1, 3, 2); hold on;
    subplot(1, 3, 3); hold on
end

% Go through series of step intensities (slow).
clear SSCur stimIntensity FlashAmp
FlashScFact = 1;
nStepIntensities = 11;
for step = 1:nStepIntensities
    
    % Create stimulus: step along.
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

    % Compute adapted current via old code.
    params.bgVolts  = 0;
    [~, adaptedCurStepOnly] = coneAdapt(sensor, 'linear', params);
    
    % Compute using outersegment object osLinear.
    noiseFlag = 0;
    adaptedOSStepOnly = osLinear();
    adaptedOSStepOnly = osSet(adaptedOSStepOnly, 'noiseFlag', noiseFlag);
    adaptedOSStepOnly = osCompute(adaptedOSStepOnly, sensor);

    % Create stimulus: step + flash.
    stimulus = zeros(nSamples, 1);
    stimulus(stimPeriod(1):stimPeriod(2)) = stimIntensity(step);
    stimulus(flashTime(1):flashTime(1)+flashDur) = stimulus(flashTime(1):flashTime(1)+flashDur) + flashIntens;
    stimulus(flashTime(2):flashTime(2)+flashDur) = stimulus(flashTime(2):flashTime(2)+flashDur) + flashIntens*FlashScFact;
    stimulus = reshape(stimulus, [1 1 nSamples]);

    % Set photon rates.
    sensor = sensorSet(sensor, 'photon rate', stimulus);

    % Compute adapted current via old code.
    params.bgVolts  = 0;
    [~, adaptedCur] = coneAdapt(sensor, 'linear', params);
     
    % Compute using outersegment object.
    noiseFlag = 0;
    adaptedOS = osLinear();
    adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag);   
    adaptedOS = osCompute(adaptedOS, sensor);

    % Plot.
    if (runTimeParams.generatePlots)
        subplot(1, 3, 1);
        plot((1:nSamples)*timeStep, adaptedCur(:), 'r', 'lineWidth', 4);
        plot((1:nSamples)*timeStep, adaptedOS.coneCurrentSignal(:), ':k', 'lineWidth', 2);

        subplot(1, 3, 2);
        plot((1:nSamples)*timeStep, adaptedCurStepOnly(:), 'r', 'lineWidth', 4);
        plot((1:nSamples)*timeStep, adaptedOSStepOnly.coneCurrentSignal(:), ':k', 'lineWidth', 2);

        subplot(1, 3, 3);
        plot((1:nSamples)*timeStep, adaptedCur(:) - adaptedCurStepOnly(:), 'r', 'lineWidth', 2);       
        plot((1:nSamples)*timeStep, adaptedOS.coneCurrentSignal(:) - adaptedOSStepOnly.coneCurrentSignal(:), ':k', 'lineWidth', 2);
    end

    % Summary measures.
    tempVar = adaptedCur(:) - adaptedCurStepOnly(:);
    tempOS = adaptedOS.coneCurrentSignal(:) - adaptedOSStepOnly.coneCurrentSignal(:); 
    flashAmp(step) = max(tempVar(flashTime(2):flashTime(2)+1000)) / (FlashScFact * max(tempVar(flashTime(1):flashTime(1)+1000)));
    ssCur(step) = -(adaptedCurStepOnly(1, 1, stimPeriod(2)) - adaptedCurStepOnly(1, 1, 1))/adaptedCurStepOnly(1, 1, 1);
    [maxVal, maxLoc] = max(tempVar(flashTime(2):flashTime(2)+1000));
    tPeak(step) = maxLoc*timeStep;
    flashAmpOS(step) = max(tempOS(flashTime(2):flashTime(2)+1000)) / (FlashScFact * max(tempOS(flashTime(1):flashTime(1)+1000)));
    ssCurOS(step) = -(adaptedOSStepOnly.coneCurrentSignal(1, 1, stimPeriod(2)) - adaptedOSStepOnly.coneCurrentSignal(1, 1, 1))/adaptedOSStepOnly.coneCurrentSignal(1, 1, 1);
    
    % Get transient of step response, computed both ways.
    transientTime = stimPeriod(1)+round(6e-2/timeStep):stimPeriod(2);
    transient{step} = adaptedCurStepOnly(transientTime);
    transient{step} = transient{step} - transient{step}(end);
    transientOS{step} = adaptedOSStepOnly.coneCurrentSignal(transientTime);
    transientOS{step} = transientOS{step} - transientOS{step}(end);
end

% Label plot.
if (runTimeParams.generatePlots)
    subplot(1, 3, 1);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('Adapted Current');
    legend('original code', 'osLinear object','Location','NorthEast');
    
    subplot(1, 3, 2);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('Adapted Current Step Only');
    legend('original code', 'osLinear object','Location','NorthEast');
    
    subplot(1, 3, 3);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('Adaoted Current Steady State');    
    legend('original code', 'osLinear object','Location','NorthEast');
    
end

% Plot sensitivity vs background intensity relation
%
% Half desens around 2500 (Angueyra and Rieke, 2013)
if (runTimeParams.generatePlots)
    vcNewGraphWin; 
    loglog(stimIntensity, flashAmp, 'ro'); 
    hold on;
    loglog(stimIntensity, flashAmpOS, 'kx');
    xlabel('Background','fontsize',16);
    ylabel('Sensitivity','fontsize',16);
    axis tight;    
    title('Sensitivity vs. Background Intensity','fontsize',16);
    legend('coneAdapt','osLinear');
end

% Plot transient from last step simulated
if (runTimeParams.generatePlots)
    vcNewGraphWin; hold on
    plot(transientTime, transient{end}(:),'r','Linewidth',4);%, tme, fit);
    plot(transientTime, transientOS{end}(:),':k','Linewidth',2);%, tme, fit);
    xlabel('Time (sec)','fontsize',16);
    ylabel('Photocurrent (pA)','fontsize',16);
    title('Transient for largest step','fontsize',16);
    legend('coneAdapt','osLinear');
end

% Save validation data
tolerance = 1e-12;
UnitTest.assertIsZero(max(abs(flashAmp(:)-flashAmpOS(:))),'Comparison for dark impulse response',tolerance);
UnitTest.validationData('adaptedCur',adaptedCur);
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
stimulus = eyeMovementExample.data.Stim;
nSamples = length(stimulus);
stimulus = reshape(stimulus, [1 1 nSamples]);

% Create human sensor.
sensor = sensorCreate('human');
sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
sensor = sensorSet(sensor, 'time interval', 5e-5);

% Set photon rates.
% This is an artifact of directly specifying the stimulus
% in the sensor, and will not be an issue when the sensor
% is the result of a sensorCompute command on a scene and oi.
sensor = sensorSet(sensor, 'photon rate', stimulus);

% Compute adapted current with coneAdapt function
params.bgVolts  = 0;
[~, adaptedCur] = coneAdapt(sensor, 'linear', params);

% Create outersegment object.
noiseFlag = 0;
adaptedOS = osLinear();
adaptedOS = osSet(adaptedOS, 'noiseFlag', noiseFlag);
adaptedOS = osCompute(adaptedOS, sensor);
osAdaptedCur = osGet(adaptedOS, 'coneCurrentSignal');

% Plot the two calculations and compare against measured data.
if (runTimeParams.generatePlots)
    vcNewGraphWin; hold on;
    % Correct for DC offset to measuredCur by setting last value to zero.
    plot((1:nSamples)*5e-5, measuredCur(:)-measuredCur(end), 'b', 'LineWidth', 4);
    plot((1:nSamples)*5e-5, adaptedCur(:)-adaptedCur(end), 'r', 'LineWidth', 3);
    plot((1:nSamples)*5e-5, osAdaptedCur(:)-osAdaptedCur(end), 'k:', 'LineWidth', 2);
    axis tight
    xlabel('Time (sec)','FontSize',14);
    ylabel('Photocurrent (pA)','FontSize',14);
    title('Comparison for end subtracted eye movement example data');
    legend('coneAdapt','osLinear');
end

% Save validation data
tolerance = 1e-11;
UnitTest.assertIsZero(max(abs(adaptedCur-osAdaptedCur)),'Comparison for dark impulse response',tolerance);
UnitTest.validationData('adaptedCur',adaptedCur);
UnitTest.validationData('osAdaptedCur',osAdaptedCur);
clear adaptedOS paramsOS

%% Increment/decrement asymmetry

% Set up parameters for stimulus.
stimPeriod = [5000 7500]; % start and end time of step
nSamples = stimPeriod(2)+5000;
timeStep = 2e-4;
contrast = 1;

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
clear maxInc maxDec stimIntensity
stepLevels = 7;
for step = 1:stepLevels
    % Create stimulus: step alone.
    stimIntensity(step) = 250 * 2^step;          % photons per second
    stimulusInc = zeros(nSamples, 1);
    stimulusInc(100:nSamples-100) = stimIntensity(step);
    stimulusInc(stimPeriod(1):stimPeriod(2)) = stimIntensity(step)*(1 + contrast);
    stimulus = reshape(stimulusInc, [1 1 nSamples]);

    % Set photon rates.
    % This is an artifact of directly specifying the stimulus
    % in the sensor, and will not be an issue when the sensor
    % is the result of a sensorCompute command on a scene and oi.
    sensor = sensorSet(sensor, 'photon rate', stimulus);

    % Compute adapted current using coneAdapt.
    params.bgVolts  = 0;
    [~, adaptedCurInc] = coneAdapt(sensor, 'linear', params);
    
    % Create outersegment object and compute with it.
    noiseFlag = 0;
    adaptedOSInc = osLinear();
    adaptedOSInc = osSet(adaptedOSInc, 'noiseFlag', noiseFlag);  
    adaptedOSInc = osCompute(adaptedOSInc, sensor);
    
    % Create stimulus: step + flashes.
    stimulusDec = zeros(nSamples, 1);
    stimulusDec(100:nSamples-100) = stimIntensity(step);
    stimulusDec(stimPeriod(1):stimPeriod(2)) = stimIntensity(step)*(1 - contrast);
    stimulus = reshape(stimulusDec, [1 1 nSamples]);

    % Set photon rates.
    sensor = sensorSet(sensor, 'photon rate', stimulus);

    % Compute adapted current with cone adapt.
    params.bgVolts  = 0;
    [~, adaptedCurDec] = coneAdapt(sensor, 'linear', params);
    
    % And with outersegment object.
    noiseFlag = 0;
    adaptedOSDec = osLinear();
    adaptedOSDec = osSet(adaptedOSDec, 'noiseFlag', noiseFlag);  
    adaptedOSDec = osCompute(adaptedOSDec, sensor);

    % Summary measures.
    maxInc(step) = adaptedCurInc(1, 1, stimPeriod(2)-1) - adaptedCurInc(1, 1, stimPeriod(1)-1);
    maxDec(step) = adaptedCurDec(1, 1, stimPeriod(2)-1) - adaptedCurDec(1, 1, stimPeriod(1)-1);
    
    maxOSInc(step) = adaptedOSInc.coneCurrentSignal(1, 1, stimPeriod(2)-1) - adaptedOSInc.coneCurrentSignal(1, 1, stimPeriod(1)-1);
    maxOSDec(step) = adaptedOSDec.coneCurrentSignal(1, 1, stimPeriod(2)-1) - adaptedOSDec.coneCurrentSignal(1, 1, stimPeriod(1)-1);
    
    % Add current to the figure.
    if (runTimeParams.generatePlots)
        subplot(1,2,1); hold on
        plot((1:nSamples)*timeStep, adaptedCurInc(:), (1:nSamples)*timeStep, adaptedCurDec(:));
        
        subplot(1,2,2); hold on
        plot((1:nSamples)*timeStep, adaptedOSInc.coneCurrentSignal(:), (1:nSamples)*timeStep, adaptedOSDec.coneCurrentSignal(:));
    end
end

% Finish off the figure
if (runTimeParams.generatePlots)
    subplot(1,2,1);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('coneAdapt');
    subplot(1,2,2);
    xlabel('Time (sec)');
    ylabel('Photocurrent (pA)');
    title('osLinear');
end

% Plot summary measure computed the two ways
if (runTimeParams.generatePlots)
    vcNewGraphWin; 
    semilogx(stimIntensity, -maxDec ./ maxInc, 'ro');
    hold on;
    semilogx(stimIntensity, -maxOSDec ./ maxOSInc, 'kx');
    xlabel('Background','fontsize',16);
    ylabel('Dec/Inc Response Ratio','fontsize',16);
    legend('coneAdapt','osLinear','location','northwest');
    title('Dec/Inc Response Ratio vs. Background Intensity','fontsize',14);
end

% Tuck data away.
tolerance = 1e-12;
UnitTest.assertIsZero(max(abs((-maxDec ./ maxInc)-(-maxOSDec ./ maxOSInc))),'Comparison for dark impulse response',tolerance);
UnitTest.validationData('adaptedCurInc',adaptedCurInc);
UnitTest.validationData('osAdaptedCurInc',adaptedOSInc.coneCurrentSignal);
UnitTest.validationData('adaptedCurDec',adaptedCurDec);
UnitTest.validationData('osAdaptedCurDec',adaptedOSDec.coneCurrentSignal);

end

% ***********************
% Original coneAdapt and riekeLinearCone functions follow.  These are no
% longer in the isetbio distribution.
% ***********************

function [sensor, adaptedData] = coneAdapt(sensor, typeAdapt, params)
% Cone adaptation
%
%   [sensor, adaptedData] = coneAdapt(sensor, typeAdapt, varargin)
%
% Implement adaptation models to produce the cone volts. Cone absorption
% sampels should be computed and stored in sensorGet(sensor, volts).
%
% The mean level is determined as a reference so that increments and
% decrements are meaningful. This is needed so that on- and off-cells can
% work. When we have a cone absorption that is above the mean, this should
% cause an off-cell to reduce its firing rate and on-cell to increase.
% Conversely, cone absorptions below the mean cause an on-cell to decrease
% and an off-cell to increase.
%
% The discussion in Rodieck book - The first steps in seeing - p. 177 and
% thereabouts for ideas about cone signals. The physiological question is
% what is the meaning of the offset or zero mean. Rodieck describes the
% effect of light as setting the mean transmitter release. In the dark,
% there is a relatively large dynamic range. As the light is on steadily,
% the release decreases and the range available for another flash. If you
% darken from the mean, the rate can increase.
%
% In addition to Rodieck's discussion, there are famous papers by Boynton
% (e.g. adaptation in cones) expressing such a model based on cone ERPs.
%
% The other issue is the total gain.  The cones can only modulate over a
% dynamic range around the mean.  So, we set the gain as well to keep the
% modulation within some range, such as +/- 40 mV of the mean (deplorizing
% voltage).
%
% The third issue is whether all the cones are the same, or there is
% space-variant adaptation. We now here consider the difference between
% different cone types, but we think it's spatial invariant.
%
% Inputs:
%  sensor     - ISETBio sensor with cone absorption computed and stored
%  typeAdapt  - The adaptation model, meaning of value are as below
%           0 = no adaptation
%           1 = a single gain map for the whole cone array
%           2 = one gain map computed for each cone class
%           3 = non-linear adaptation
%           4 = cone adaptation by physiological differential equation
%               (Rieke)
%  params     - Cone adaptation parameters, could include
%    .bgVolts  = background voltage
%    .vSwing   = maximum volts response
%    .addNoise = boolean, used in the Rieke cone adapt case
%
%
% Output:
%   sensor      - ISETBio sensor with cone adaptation parameters set. The
%                 parameters include gain and offset. You could retrieve
%                 the adaptation parameters and data by sensorGet function
%                 calls.
%   adaptedData - Adapted voltage 3D matrix, would be in same size as volts
%                 image in sensor.
%
% Notes:
%   In the future, we will look at the time series and have a time-varying
%   adaptation function. Generally, cones take 200ms and rods take 800ms
%   for adaptation.
%
% Examples:
%   scene = sceneCreate; oi = oiCreate; oi = oiCompute(oi,scene);
%   sensor = sensorCreate('human'); sensor = sensorCompute(sensor,oi);
%
%   sensor = coneAdapt(sensor, 1);
%   adaptedData = sensorGet(sensor, 'adapted volts');
%   gain = sensorGet(sensor, 'adaptation gain');
%   offset = sensorGet(sensor, 'adaptation volts');
%
%   sensor = coneAdapt(sensor,4);  % Rieke model 
%   v = sensorGet(sensor,'volts');
%
% (c) Stanford VISTA Lab, 2014

%% Warn of imminent squashing
% warning('coneAdapt will soon be deprecated and replaced by the outersegment class; see t_OuterSegmentClasses.m, v_osLinearObjectTemp.m and v_osBioPhysObjectTemp.m');

%% Check inputs and Init
if notDefined('sensor'),      error('sensor is required'); end
if ~sensorCheckHuman(sensor), error('unsupported species'); end
if notDefined('typeAdapt'),   typeAdapt = 3; end
if notDefined('params'),      params = []; end

%% Compute cone adaptations
volts  = double(sensorGet(sensor, 'volts'));
if isempty(volts), error('cone absorptions should be pre-computed'); end

if isfield(params, 'vSwing'),    vSwing = params.vSwing;
else                             vSwing = sensorGet(sensor,'pixel voltageSwing');
end

if isfield(params, 'bgVolts'),     bgVolts = params.bgVolts;
else                               bgVolts = mean(volts(:));
end

if isfield(params, 'addNoise'),  noiseFlag = params.addNoise;
else                             noiseFlag = false;
end
    

if ischar(typeAdapt), typeAdapt = ieParamFormat(typeAdapt); end
    
switch typeAdapt
    case {0, 'noadapt', 'noadaptation'} % no adaptation
        gainMap = 1;
        bgVolts  = 0;
        
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);
        
    case {1, 'globalgain'}
        % Use same gain for all cone type
        
        % Adjust for bg volts as offset
        volts = volts - bgVolts;
        
        % Set the gain so that the max - min get scaled to vSwig mV
        gainMap = vSwing / (max(volts(:)) - min(volts(:)));
        adaptedData = gainMap * volts;
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);
        
    case {2, 'adaptbyconetype'}
        % Use different gains for each cone type
        % For human, the cone types in sensor are 1~4 for K, L, M, S
        if isscalar(bgVolts), bgVolts = ones(4,1)*bgVolts; end
        assert(numel(bgVolts) == 4, 'bgVolts should be of 4 elements');
        gainMap = ones(4, 1);
        
        % Adjust for backgroud
        volts = volts - bgVolts(sensorGet(sensor, 'cone type'));
        
        % Compute gain map
        for ii = 2 : 4 % L,M,S and we don't need to compute for K
            v = sensorGet(sensor,'volts',ii);
            if ~isempty(v)
                gainMap(ii) = vSwing / (max(v) - min(v));
            end
        end
        
        nSamples = size(volts, 3);
        gainMap = gainMap(sensorGet(sensor, 'cone type'));
        
        adaptedData = volts .* repmat(gainMap, [1 1 nSamples]);
        
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);

    case {3, 'felice'}
        % In this case, we will do non-linear cone adaptation
        % Reference:
        %   Felice A. Dunn, et al. Light adaptation in cone vision involves
        %   switching between receptor and post-receptor sites,
        %   doi:10.1038/nature06150
        % From there and some other papers, we know two facts:
        %   1) The steady state response follows:
        %      R/R_max = 1/(1 + (I_{1/2}/I_B)^n)
        %   2) The cone adaptation follows
        %      AMP/AMP_{dark} = ((a + bI_B)/(a + I_B))/(cI_B+1)
        % To compute the adapted voltage for input I on background I_B, we
        % use the following process
        %   1) Compute steady state response for I as R_i
        %   2) Compute dark adapted response for I as R_i * AMP / AMP{dark}
        %   3) Compute I_B adapted response
        %
        % Notes:
        %   In the original paper, the formulas are only tested and fitted
        %   for L cones. Here, we use it for all cone types and this might
        %   lead to simulation errors.
        %
        %   Also, by this method, the adapted cone volts are still positive
        %   This is different from the previous two method
        %
        % (HJ) ISETBIO TEAM, 2014
        
        % Check inputs
        if ~isscalar(bgVolts), error('bgVolts should be scalar'); end
        
        % Init function handles
        cg = 1/pixelGet(sensorGet(sensor, 'pixel'), 'conversion gain');
        cg = cg / sensorGet(sensor, 'exp time');
        steadyR  = @(x) 1./(1+(45000 ./ x / cg).^0.7)*vSwing;
        ampRatio = @(x) (100 + 1.3*x/cg)./(100 + x/cg)./(0.00029*x/cg + 1);
        
        % Compute steady state response
        sR = steadyR(volts);
        
        % Compute dark adapted response
        dR = sR ./ ampRatio(volts);
        
        % Compute adapted response at level I_B
        adaptedData = dR * ampRatio(bgVolts);
        
        % Set gain map and offset
        gainMap = adaptedData ./ volts;
        bgVolts = 0; % We don't have an actual offset
        
        % Set adaptation parameters back to sensor
        sensor = sensorSet(sensor, 'adaptation gain', gainMap);
        sensor = sensorSet(sensor, 'adaptation offset', bgVolts);

    case {4, 'rieke'}
        % See rieke<TAB> for explanations.
        
        % sensor = sensorCreate('human');
        % sensor = sensorSet(sensor,'pixel voltage swing',0.05);
        % v = rand(32,32,200)*sensorGet(sensor,'pixel voltage swing');
        % sensor = sensorSet(sensor,'volts',v);
        % vcAddObject(sensor); sensorWindow;
        % [~,adaptedData] = coneAdapt(sensor, 4);
        % s = 255 / max(abs(adaptedData(:)));
        %                 vcNewGraphWin
        %                 for ii=1:size(adaptedData,3)
        %                     imagesc(abs(s*adaptedData(:,:,ii)));
        %                     pause(0.02);
        %                 end
        % 
        % Fix this mplay(adaptedData,'I');
        
        p = osInit;
        expTime = sensorGet(sensor,'exposure time');
        sz = sensorGet(sensor,'size');
        
        % absRate = sensorGet(sensor,'absorptions per second');        
        % pRate = sensorGet(sensor, 'volts') / (sensorGet(sensor,'conversion gain')*expTime);
        pRate = sensorGet(sensor, 'photon rate');
        
        % Compute background adaptation parameters
        bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);
        
        initialState = osAdaptSteadyState(bgR, p, sz);
        initialState.timeInterval = sensorGet(sensor, 'time interval');
        adaptedData  = osAdaptTemporal(pRate, initialState);
        
        if noiseFlag
            adaptedData = osAddNoise(adaptedData, params);
        end
  
    case {5, 'linear'}
        expTime = sensorGet(sensor,'exposure time');
        sz = sensorGet(sensor,'size');
        
        % absRate = sensorGet(sensor,'absorptions per second');        
        % pRate2 = sensorGet(sensor, 'volts') / (sensorGet(sensor,'conversion gain')*expTime);
        pRate = sensorGet(sensor, 'photon rate');  
        
        % Compute background adaptation parameters
        bgR = 0;
        
        initialState = osInit;
        initialState.timeInterval = sensorGet(sensor, 'time interval');
        if isfield(params, 'Compress')
            initialState.Compress = params.Compress;
        else
            initialState.Compress = 0;
        end
        if isfield(params, 'Ib')
            initialState.Ib = params.Ib;
        end
        adaptedData  = riekeLinearCone(pRate, initialState);
        if noiseFlag
            adaptedData = osAddNoise(adaptedData);
        end
       
    otherwise
        error('unknown adaptation type');
end



end


function adaptedData = riekeLinearCone(pRate,p)

%%
if ~exist('pRate','var') || isempty(pRate), error('Photon absorption rate required.'); end
if ~exist('p','var') || isempty(p), error('Initial state required.'); end
    
% One millisecond time step, as in s_fredConeModel.m
if isfield(p, 'timeInterval'), dt = p.timeInterval; else dt = 0.001; end

if isfield(p, 'Ib'), Ib = p.Ib; else Ib = mean(pRate(:)); end

if isfield(p, 'Compress'), CompressFlag = p.Compress; else CompressFlag = 0; end

%% cone IRF


if (ndims(pRate) == 3)
    timeBins = size(pRate, 3);
elseif (ndims(pRate) == 2)
    timeBins = size(pRate, 2);
end

% only compute filter for 0.3 seconds
TimeAxis = (1:timeBins)*dt;
TimeAxis = TimeAxis(find(TimeAxis <= 0.3));

% move these to riekeInit
TauR = 0.0216;   % Rising Phase Time Constant
TauD = 0.0299;   % Damping Time Constant
TauP = 0.5311;   % Period

Io = 2250;                     % half-desensitizing background (in R*/cone/sec, from Juan's paper - corrected)if (i
gain_dark = 0.22; %0.32;              % from Juan's paper (approximate peak of the IRF measured in darkness, and in units of pA/R*) - corrected
gainRatio = 1 ./ (1+(Ib./Io)); % the right side of the equation above, and the gain ratio implied by the bkgnd adapting field

Phi = 34.1814;   % Phase

Filter = (((TimeAxis./TauR).^3)./(1+((TimeAxis./TauR).^3))) .* exp(-((TimeAxis./TauD))).*cos(((2.*pi.*TimeAxis)./TauP)+(2*pi*Phi/360));

newGain = gainRatio .* gain_dark;
oldGain = max(Filter);
Filter = Filter * dt * newGain / oldGain;
maxCur = p.k * p.gdark^p.h/2;
meanCur = maxCur * (1 - 1 / (1 + 45000 / mean(pRate(:))));
FilterFFT = fft(Filter');

if (ndims(pRate) == 3)
    % pre-allocate memory
    adaptedData = zeros(size(pRate,1), size(pRate, 2), timeBins); 
    
    for x = 1:size(pRate, 1)
        for y = 1:size(pRate, 2)
            tempData = conv(squeeze(pRate(x, y, :)), Filter);
    %        tempData = real(ifft(conj(fft(squeeze(pRate(x, y, :))) .* FilterFFT)));
            if (CompressFlag)
                tempData = tempData / maxCur;
                tempData = meanCur * (tempData ./ (1 + 1 ./ tempData)-1);
            else
                tempData = tempData - meanCur;
            end
            adaptedData(x, y, :) = tempData(1:timeBins);
        end
    end
    
elseif (ndims(pRate) == 2)
    
    % pre-allocate memory
    adaptedData = zeros(size(pRate,1),timeBins); 
    
    for xy = 1:size(pRate, 1)
        tempData = conv(squeeze(pRate(xy, :)), Filter);
        if (p.Compress)
            tempData = tempData / maxCur;
            tempData = meanCur * (tempData ./ (1 + 1 ./ tempData)-1);
        else
            tempData = tempData - meanCur;
        end
        adaptedData(xy, :) = tempData(1:timeBins); 
    end
end

end
