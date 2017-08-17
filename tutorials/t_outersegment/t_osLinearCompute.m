function t_osLinearCompute
%
% Description:
%    Demonstrates how to compute photocurrent responses for a mosaic that
%    employs the @osLinear outer-segment model, in response to a dynamic stimulus
%    whose contrast is modulated against an adapting background. 
%
% 
% NPC, ISETBIO Team, 2017

% Basic parameters
eccDegs = 5;            % mosaic eccentricity in visual degrees (> 10 results in peripheral os dynamics)
fov = 0.5;              % field of view in degrees
meanLuminance = 50;     % stimulus mean luminance (cd/m2)
nTrials = 1;            % response instances to compute
testContrast = [0.1];   % stimulus contrasts to examine

% Set the random seed
rng(1);

% Default params for different components
[spatialParams, temporalParams, colorModulationParams, backgroundParams] = getStimParams(fov, meanLuminance);
[oiParams, mosaicParams] = getPhysioParams(fov, eccDegs);

%% Generate the cone mosaic
cm = coneMosaic('os',osLinear('eccentricity', eccDegs)); % creating the cone mosaic, currently on hold 
% Set the field of view
cm.setSizeToFOV(mosaicParams.fov);
% Set the integration time
cm.integrationTime = mosaicParams.integrationTimeInSeconds;
% Set the isomerization noise. choose from {'random',  'frozen', 'none'}
cm.noiseFlag = 'frozen';           
% Set the photocurrent noise. choose from {'random',  'frozen', 'none'}
cm.os.noiseFlag = 'frozen';

%% Create the background scene (zero contrast)
theBaseColorModulationParams = colorModulationParams;
theBaseColorModulationParams.coneContrasts = [0 0 0]';
theBaseColorModulationParams.contrast = 0;
backgroundScene = colorSceneCreate(spatialParams, backgroundParams, theBaseColorModulationParams, []);

%% Compute the background OI
BaseOI = colorDetectOpticalImageConstruct(oiParams);
oiBackground = BaseOI;
oiBackground = oiCompute(oiBackground, backgroundScene);
oiModulated = oiBackground;

%% Compute the background OI sequence
zeroContrastOIsequence = oiSequence(oiBackground, oiModulated, temporalParams.sampleTimes, temporalParams.TemporalWindow, 'composition', 'blend');

%% Generate eye movement path for zero contrast stimulus
eyeMovementsNum = zeroContrastOIsequence.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
theZeroContrastEMpaths = colorDetectMultiTrialEMPathGenerate(cm, 1, eyeMovementsNum, 'none');


%% Compute response to adapting stimulus, the outer segment IR functions and the mean photocurrents (to the adapting stimulus)
% The outer segment IR functions and mean photocurrents to the adapting
% stimulus are used in the later computation of the response to the test stimuli
[isomerizationsAdaptingStim, photocurrentsAdaptingStim, osImpulseResponseFunctions, osMeanCurrents] = ...
     cm.computeForOISequence(zeroContrastOIsequence, ...
        'emPaths', theZeroContrastEMpaths, ...
        'interpFilters', [], ...
        'meanCur', [], ...
        'currentFlag', true);

%% Compute responses for all contrasts examined
for iContrast = 1:numel(testContrast)
    %% Create the modulated scene
    colorModulationParams.contrast = testContrast(iContrast);
    modulatedScene = colorSceneCreate(spatialParams, backgroundParams, colorModulationParams, []);

    %% Compute the modulated OI
    oiModulated = BaseOI;
    oiModulated = oiCompute(oiModulated, modulatedScene);
    
    %% Compute the stimulus oi sequence
    stimulusOIsequence = oiSequence(oiBackground, oiModulated, temporalParams.sampleTimes, temporalParams.TemporalWindow, 'composition', 'blend');
    %%  Visualize the stimuus oiSequence
	%stimulusOIsequence.visualize('format', 'montage', 'showIlluminanceMap', true);
    
    %% Generate eye movement paths for all response instances
    eyeMovementsNum = stimulusOIsequence.maxEyeMovementsNumGivenIntegrationTime(cm.integrationTime);
    theEMpaths = colorDetectMultiTrialEMPathGenerate(cm, nTrials, eyeMovementsNum, temporalParams.emPathType);

    %% Compute noise-free responses - this is only used for improving the response visualization. Otherwise it can be skipped
    noiseFlags = {cm.noiseFlag, cm.os.noiseFlag};
    cm.noiseFlag = 'none';
    cm.os.noiseFlag = 'none';
    [isomerizationsNoiseFree{iContrast}, photocurrentsNoiseFree{iContrast}, ~, ~] = ...
        cm.computeForOISequence(stimulusOIsequence, ...
                    'emPaths', 0*theEMpaths(1,:,:), ...    % zero movement path for the noise-free versions
                    'interpFilters', osImpulseResponseFunctions, ...
                    'meanCur', osMeanCurrents, ...
                    'currentFlag', true);
     
    %% Compute the noisy responses
    cm.noiseFlag = noiseFlags{1};
    cm.os.noiseFlag = noiseFlags{2};
    [isomerizations{iContrast}, photocurrents{iContrast}, ~, ~] = ...
        cm.computeForOISequence(stimulusOIsequence, ...
                    'emPaths', theEMpaths, ...
                    'interpFilters', osImpulseResponseFunctions, ...
                    'meanCur', osMeanCurrents, ...
                    'currentFlag', true);
end  % iContrast 


%% Visualize responses
trialVisualized = 1;
maxConesVisualized = 5000;
contrastVisualized = 1;
photocurrentRange = [-100 -30];
meanLMPhotocurrentRange = 1.2*[-1 1];
meanSPhotocurrentRange = 1.2*[-1 1];

% Visualize photocurrents
visualizeResponses(cm, trialVisualized, maxConesVisualized, contrastVisualized, ...
    photocurrentsNoiseFree, photocurrents, photocurrentsAdaptingStim, osMeanCurrents, ...
    photocurrentRange, meanLMPhotocurrentRange, meanSPhotocurrentRange);

end


%% Visualize routine
function visualizeResponses(cm, trialVisualized, maxConesVisualized, contrastVisualized, ...
    photocurrentsNoiseFree, photocurrents, photocurrentsAdaptingStim, osMeanCurrents,...
    photocurrentRange, meanLMPhotocurrentRange, meanSPhotocurrentRange)

currNoiseFree = RGB2XWFormat(squeeze(photocurrentsNoiseFree{contrastVisualized}(trialVisualized,:,:,:)));
curr = RGB2XWFormat(squeeze(photocurrents{contrastVisualized}(trialVisualized,:,:,:)));
currAdaptingStim = RGB2XWFormat(squeeze(photocurrentsAdaptingStim(1,:,:,:)));

[peakIncLconeIndices, peakDecLconeIndices] = peakResponseConeIndices(cm, currNoiseFree, osMeanCurrents, 'L');
[peakIncMconeIndices, peakDecMconeIndices] = peakResponseConeIndices(cm, currNoiseFree, osMeanCurrents, 'M');
[peakIncSconeIndices, peakDecSconeIndices] = peakResponseConeIndices(cm, currNoiseFree, osMeanCurrents, 'S');

incLconesNumVisualized = min([numel(peakIncLconeIndices) maxConesVisualized]);
incMconesNumVisualized = min([numel(peakIncMconeIndices) maxConesVisualized]);
incSconesNumVisualized = min([numel(peakIncSconeIndices) maxConesVisualized]);

decLconesNumVisualized = min([numel(peakDecLconeIndices) maxConesVisualized]);
decMconesNumVisualized = min([numel(peakDecMconeIndices) maxConesVisualized]);
decSconesNumVisualized = min([numel(peakDecSconeIndices) maxConesVisualized]);

peakIncLconeResponsesNoiseFree = currNoiseFree(peakIncLconeIndices(1:incLconesNumVisualized),:);
peakIncMconeResponsesNoiseFree = currNoiseFree(peakIncMconeIndices(1:incMconesNumVisualized),:);
peakIncSconeResponsesNoiseFree = currNoiseFree(peakIncSconeIndices(1:incSconesNumVisualized),:);

peakDecLconeResponsesNoiseFree = currNoiseFree(peakDecLconeIndices(1:decLconesNumVisualized),:);
peakDecMconeResponsesNoiseFree = currNoiseFree(peakDecMconeIndices(1:decMconesNumVisualized),:);
peakDecSconeResponsesNoiseFree = currNoiseFree(peakDecSconeIndices(1:decSconesNumVisualized),:);

peakIncLconeResponses = curr(peakIncLconeIndices(1:incLconesNumVisualized),:);
peakIncMconeResponses = curr(peakIncMconeIndices(1:incMconesNumVisualized),:);
peakIncSconeResponses = curr(peakIncSconeIndices(1:incSconesNumVisualized),:);

peakDecLconeResponses = curr(peakDecLconeIndices(1:decLconesNumVisualized),:);
peakDecMconeResponses = curr(peakDecMconeIndices(1:decMconesNumVisualized),:);
peakDecSconeResponses = curr(peakDecSconeIndices(1:decSconesNumVisualized),:);

peakIncLconeResponsesAdaptingStim = currAdaptingStim(peakIncLconeIndices(1:incLconesNumVisualized),:);
peakIncMconeResponsesAdaptingStim = currAdaptingStim(peakIncMconeIndices(1:incMconesNumVisualized),:);
peakIncSconeResponsesAdaptingStim = currAdaptingStim(peakIncSconeIndices(1:incSconesNumVisualized),:);

peakDecLconeResponsesAdaptingStim = currAdaptingStim(peakDecLconeIndices(1:decLconesNumVisualized),:);
peakDecMconeResponsesAdaptingStim = currAdaptingStim(peakDecMconeIndices(1:decMconesNumVisualized),:);
peakDecSconeResponsesAdaptingStim = currAdaptingStim(peakDecSconeIndices(1:decSconesNumVisualized),:);

timeAxis = cm.timeAxis;

LincColor = [1 0.3 0.3];
LdecColor = [0.6 0 0];

MincColor = [0.3 0.8 0.3];
MdecColor = [0.0 0.4 0.0];

SincColor = [0.3 0.3 1];
SdecColor = [0.0 0 0.6];

subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 6, ...
       'heightMargin',  0.08, ...
       'widthMargin',    0.03, ...
       'leftMargin',     0.04, ...
       'rightMargin',    0.0, ...
       'bottomMargin',   0.06, ...
       'topMargin',      0.03);
   
hFig = figure();
clf;
set(hFig, 'Position', [10 10 1600 700], 'Color', [1 1 1]);
subplot('Position', subplotPosVectors(1,1).v);
plot(timeAxis, peakIncLconeResponsesAdaptingStim', 'r-'); hold on;
plot(timeAxis, peakDecLconeResponsesAdaptingStim', 'r-');
set(gca, 'YLim', photocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:5:0, 'FontSize', 14);
grid on; box off;
ylabel('photocurrent (pAmps)');
title('L cone responses (adapt)');

subplot('Position', subplotPosVectors(1,2).v);
plot(timeAxis, peakIncLconeResponses', 'r-'); hold on;
plot(timeAxis, peakDecLconeResponses', 'r-'); 
set(gca, 'YLim', photocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:5:0, 'FontSize', 14);
grid on; box off;
title('L cone responses (test)');

subplot('Position', subplotPosVectors(1,3).v);
plot(timeAxis, peakIncMconeResponsesAdaptingStim', '-', 'Color', [0 0.8 0.0]); hold on;
plot(timeAxis, peakDecMconeResponsesAdaptingStim', '-', 'Color', [0 0.8 0.0]);
set(gca, 'YLim', photocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:5:0, 'FontSize', 14);
grid on; box off;
title('M cone responses (adapt)');

subplot('Position', subplotPosVectors(1,4).v);
plot(timeAxis, peakIncMconeResponses', '-', 'Color', [0 0.7 0.0]); hold on;
plot(timeAxis, peakDecMconeResponses', '-', 'Color', [0 0.7 0.0]);
set(gca, 'YLim', photocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:5:0, 'FontSize', 14);
grid on; box off;
title('M cone responses (test)');

subplot('Position', subplotPosVectors(1,5).v);
plot(timeAxis, peakIncSconeResponsesAdaptingStim', 'b-'); hold on;
plot(timeAxis, peakDecSconeResponsesAdaptingStim', 'b-');
set(gca, 'YLim', photocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:5:0, 'FontSize', 14);
grid on; box off;
title('S cone responses (adapt)');

subplot('Position', subplotPosVectors(1,6).v);
plot(timeAxis, peakIncSconeResponses', 'b-'); hold on;
plot(timeAxis, peakDecSconeResponses', 'b-'); 
set(gca, 'YLim', photocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:5:0, 'FontSize', 14);
grid on; box off;
title('S cone responses (test)');

subplot('Position', subplotPosVectors(2,1).v);
plot(timeAxis, mean(peakIncLconeResponsesAdaptingStim,1), '-', 'Color', LincColor, 'LineWidth', 1.5); hold on;
plot(timeAxis, mean(peakDecLconeResponsesAdaptingStim,1), '-', 'Color', LdecColor, 'LineWidth', 1.5);
set(gca, 'YLim', osMeanCurrents(1) + meanLMPhotocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:1:0, 'FontSize', 14);
grid on; box off;
ylabel('photocurrent (pAmps)');
xlabel('time (seconds)');
title('mean L cone response (adapt)');

subplot('Position', subplotPosVectors(2,2).v);
plot(timeAxis, mean(peakIncLconeResponses,1), '-', 'Color', LincColor, 'LineWidth', 1.5); hold on;
plot(timeAxis, mean(peakDecLconeResponses,1), '-', 'Color', LdecColor, 'LineWidth', 1.5);
plot(timeAxis, mean(peakIncLconeResponsesNoiseFree,1), '--', 'Color', LincColor, 'LineWidth', 1.5);
plot(timeAxis, mean(peakDecLconeResponsesNoiseFree,1), '--', 'Color', LdecColor, 'LineWidth', 1.5);
set(gca, 'YLim', osMeanCurrents(1) + meanLMPhotocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:1:0, 'FontSize', 14);
legend({...
    sprintf('mean-inc (%d)', incLconesNumVisualized), ...
    sprintf('mean-dec (%d)', decLconesNumVisualized), ...
    'mean-inc (noise-free)', 'mean-dec (noise-free)'}, 'Location', 'South');
grid on; box off;
title('mean L cone response (test)');

subplot('Position', subplotPosVectors(2,3).v);
plot(timeAxis, mean(peakIncMconeResponsesAdaptingStim,1), '-', 'Color', MincColor, 'LineWidth', 1.5); hold on;
plot(timeAxis, mean(peakDecMconeResponsesAdaptingStim,1), '-', 'Color', MdecColor, 'LineWidth', 1.5);
set(gca, 'YLim', osMeanCurrents(2) + meanLMPhotocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:1:0, 'FontSize', 14);
grid on; box off;
title('mean M cone response (adapt)');

subplot('Position', subplotPosVectors(2,4).v);
plot(timeAxis, mean(peakIncMconeResponses,1), '-', 'Color', MincColor, 'LineWidth', 1.5); hold on;
plot(timeAxis, mean(peakDecMconeResponses,1), '-', 'Color', MdecColor, 'LineWidth', 1.5);
plot(timeAxis, mean(peakIncMconeResponsesNoiseFree,1), '--', 'Color', MincColor, 'LineWidth', 1.5);
plot(timeAxis, mean(peakDecMconeResponsesNoiseFree,1), '--', 'Color', MdecColor, 'LineWidth', 1.5)
set(gca, 'YLim', osMeanCurrents(2) + meanLMPhotocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:1:0, 'FontSize', 14);
legend({...
    sprintf('mean-inc (%d)', incMconesNumVisualized), ...
    sprintf('mean-dec (%d)', decMconesNumVisualized), ...
    'mean-inc (noise-free)', 'mean-dec (noise-free)'}, 'Location', 'South');

grid on; box off;
title('mean M cone response (test)');

subplot('Position', subplotPosVectors(2,5).v);
plot(timeAxis, mean(peakIncSconeResponsesAdaptingStim,1), '-', 'Color', SincColor, 'LineWidth', 1.5); hold on;
plot(timeAxis, mean(peakDecSconeResponsesAdaptingStim,1), '-', 'Color', SdecColor, 'LineWidth', 1.5);
set(gca, 'YLim', osMeanCurrents(3) + meanSPhotocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:1:0, 'FontSize', 14);
grid on; box off;
title('mean S cone response (adapt)');

subplot('Position', subplotPosVectors(2,6).v);
plot(timeAxis, mean(peakIncSconeResponses,1), '-', 'Color', SincColor, 'LineWidth', 1.5); hold on;
plot(timeAxis, mean(peakDecSconeResponses,1), '-', 'Color', SdecColor, 'LineWidth', 1.5);
plot(timeAxis, mean(peakIncSconeResponsesNoiseFree,1), '--', 'Color', SincColor, 'LineWidth', 1.5);
plot(timeAxis, mean(peakDecSconeResponsesNoiseFree,1), '--', 'Color', SdecColor, 'LineWidth', 1.5);
set(gca, 'YLim', osMeanCurrents(3) + meanSPhotocurrentRange, 'XLim', [timeAxis(1) timeAxis(end)], 'XTick', 0:0.1:1.0, 'YTick', -100:1:0, 'FontSize', 14);
legend({...
    sprintf('mean-inc (%d)', incSconesNumVisualized), ...
    sprintf('mean-dec (%d)', decSconesNumVisualized), ...
    'mean-inc (noise-free)', 'mean-dec (noise-free)'}, 'Location', 'South');
grid on; box off;
title('mean S cone response (test)');
end        

function [peakIncrementConeIndices, peakDecrementConeIndices] = peakResponseConeIndices(cm, currents, osMeanCurrents, coneType)
% Determine target cone index
[~, coneIndex] = ismember(coneType, {'L', 'M', 'S'});
if coneIndex == 0
    peakIncrementConeIndices = [];
    peakDecrementConeIndices = [];
    return;
end

% Compute response modulations for the target cone type
pattern = cm.pattern(:);
targetConeIndices = find(pattern==coneIndex+1);
currentModulations = currents(targetConeIndices,:) - osMeanCurrents(coneIndex);
thresholdModulation = 0.1*max(abs(currentModulations(:)));

% Sort cone with positive modulation
m = max(currentModulations, [], 2);
[m, idx] = sort(m, 'descend');
idx = idx(m>thresholdModulation);
peakIncrementConeIndices = targetConeIndices(idx);

% Sort cone with negative modulation
m = min(currentModulations, [], 2);
[m, idx] = sort(m, 'ascend');
idx = idx(m<-thresholdModulation);
peakDecrementConeIndices = targetConeIndices(idx);
end


function [sampleTimes, TemporalWindow] = flatTopGaussianWindowCreate(temporalParams)
sampleTimes = 0:temporalParams.stimulusSamplingIntervalInSeconds:temporalParams.stimulusDurationInSeconds;
L = numel(sampleTimes);
upRampTime = temporalParams.windowTauInSeconds;
upRampSamples = ceil(upRampTime/temporalParams.stimulusSamplingIntervalInSeconds);
downRampSamples = L-upRampSamples+1;
TemporalWindow = ones(L,1);
TemporalWindow(1:upRampSamples) = 1-exp(-0.5*(linspace(0,1,upRampSamples)/0.33).^2);
TemporalWindow(downRampSamples:L) = TemporalWindow(upRampSamples:-1:1);
visualizeWindow = false;
if (visualizeWindow)
    figure(2);
    plot(sampleTimes, TemporalWindow, 'ks-');
    drawnow;
end
end


function [spatialParams, temporalParams, colorModulationParams,  backgroundParams] = getStimParams(fov, meanLuminance)
% getting the gabor parameters, see spatialParamsGenerate.m
spatialParams.type = 'Spatial';
spatialParams.spatialType = 'Gabor';
spatialParams.windowType = 'Gaussian';
spatialParams.fieldOfViewDegs = fov;
spatialParams.gaussianFWHMDegs = 0.942/2; % 2.3548 x standard deviation = 0.4 
spatialParams.cyclesPerDegree = 1.5*2;
spatialParams.row = 128;
spatialParams.col = 128;
spatialParams.ang = 0;
spatialParams.ph = 0;
spatialParams.viewingDistance = 0.7; % in metres, check from t_colorGaborConeAbsorptionMovie.m 

% getting the color modulation parameters, see colorModulationParamsGenerate.m 
colorModulationParams.type = 'ColorModulation';
colorModulationParams.modulationType = 'monitor';
colorModulationParams.contrast = 1;
colorModulationParams.coneContrasts = [1 1 1]'; 
colorModulationParams.coneContrasts = colorModulationParams.coneContrasts./norm(colorModulationParams.coneContrasts); % normlizing the contrast
colorModulationParams.startWl = 380;
colorModulationParams.endWl = 780;
colorModulationParams.deltaWl = 4;

% getting the background parameters, see backgroundParamsGenerate.m
lumFactor = 5.0;
backgroundParams.type = 'Background';
backgroundParams.backgroundType = 'monitor';
backgroundParams.backgroundxyY = [0.27 0.30 meanLuminance/lumFactor]'; 
backgroundParams.monitorFile = 'CRT-MODEL';
backgroundParams.leakageLum = 1.0;
backgroundParams.lumFactor = lumFactor;

% Obtaining the temporal parameters, see temporalParamsGenerate.m
temporalParams.type = 'Temporal';
temporalParams.frameRate = 60;
temporalParams.windowTauInSeconds = 0.1;
temporalParams.stimulusDurationInSeconds = 0.5; % 1000 ms stimulus presentation
temporalParams.stimulusSamplingIntervalInSeconds = 1/temporalParams.frameRate;
temporalParams.emPathType = 'none';  % choose from {'random',  'frozen', 'none'}
[temporalParams.sampleTimes,temporalParams.TemporalWindow] = flatTopGaussianWindowCreate(temporalParams); 
temporalParams.nSampleTimes = length(temporalParams.sampleTimes);
end

function [oiParams, mosaicParams] = getPhysioParams(fov, eccDegs)
% Obtaining optical image params, see oiParamsGenerate.m
oiParams.type = 'Optics';
oiParams.offAxis = false;
oiParams.blur = false; % no blurring 
oiParams.lens = false; % no effect of lens transmittance
oiParams.pupilDiamMm = 3;
oiParams.opticsModel = 'WvfHuman';
oiParams.fieldOfViewDegs = fov;

% Obtaining the cone mosaic params, see mosaicParamsGenerate.m
mosaicParams.type = 'Mosaic';
mosaicParams.conePacking = 'rect';
mosaicParams.fov = fov;
mosaicParams.realisticSconeSubmosaic = false;       % if this is set to true, there will be a 0.3 deg S-cone free region and the S-cone lattice will be semiregular
mosaicParams.LMSRatio = [0.62 0.31 0.07];
mosaicParams.innerSegmentSizeMicrons = 1.4;       
mosaicParams.apertureBlur = false;                  
mosaicParams.coneSpacingMicrons = 2.0;             
mosaicParams.mosaicRotationDegs = 0;
mosaicParams.macular = false;
mosaicParams.eccentricityDegs = eccDegs;
mosaicParams.integrationTimeInSeconds = 5/1000;
mosaicParams.osTimeStepInSeconds = 0.1/1000;
end
