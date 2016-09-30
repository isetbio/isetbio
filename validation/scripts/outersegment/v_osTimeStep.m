function varargout = v_osTimeStep(varargin)
%
% Validate current computation using different time steps. 
%
% This script tests ....
%
% 9/29/2016    npc   Created.

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;

% Generate optics
noOptics = false;
if (noOptics)
    theOI = oiCreate('diffraction limited');
    optics = oiGet(theOI,'optics');           
    optics = opticsSet(optics,'fnumber',0);
    optics = opticsSet(optics, 'off axis method', 'skip');
    theOI = oiSet(theOI,'optics', optics);
else
    theOI = oiCreate('human');
end

% Generate a uniform field scene with desired mean luminance
FOV = 2.0;  % deg
meanLuminance = 500;
theScene = uniformFieldSceneCreate(FOV, meanLuminance);
meanSceneLuminance = sceneGet(theScene, 'mean luminance');
        
%% Compute the optical image
theOI = oiCompute(theOI, theScene);
oiWavelengthAxis = oiGet(theOI, 'wave');

%% Set the stimulusSamplingInterval equal to the refresh cycle of a 120Hz display
stimulusSamplingInterval = 1/60;
stimulusTimeAxis =  -0.6:stimulusSamplingInterval:0.3;

%% The stimulus luminance is ramped up and down in time using a Gaussian window.
stimulusRampTau = 0.165;
stimulusRamp = exp(-0.5*(stimulusTimeAxis/stimulusRampTau).^2);

%% Instead of computing a new scene for each stimulus frame, we modulate the 
%% photons field of the optical image
stimulusContrast = 0.1;
backgroundPhotons = oiGet(theOI, 'photons');
retinalPhotonsTimeCourse = zeros([size(backgroundPhotons) numel(stimulusTimeAxis)]);
refRow = round(size(retinalPhotonsTimeCourse,1)/2);
refCol = round(size(retinalPhotonsTimeCourse,2)/2);
referencePositionPhotons = zeros(numel(oiWavelengthAxis), numel(stimulusTimeAxis));
for stimFrameIndex = 1:numel(stimulusTimeAxis)
    retinalPhotonsAtCurrentFrame = backgroundPhotons * (1.0 + stimulusContrast*stimulusRamp(stimFrameIndex));
    theOIsequence{stimFrameIndex} = oiSet(theOI, 'photons', retinalPhotonsAtCurrentFrame);
    referencePositionPhotons(:, stimFrameIndex) = squeeze(retinalPhotonsAtCurrentFrame(refRow, refCol, :));
end

%% Generate a human cone mosaic with 1L, 1M and 1S cone
theConeMosaic = coneMosaic;
theConeMosaic.rows = 1;
theConeMosaic.cols = 3;
theConeMosaic.pattern = [2 3 4];
theConeMosaic.noiseFlag = false;
theConeMosaic.integrationTime = stimulusSamplingInterval;


%% Generate the outer-segment object to be used by the coneMosaic
theOuterSegment = osLinear();
theOuterSegment.noiseFlag = false;
theOuterSegment.timeStep  = 1/1000;

%% Couple the outersegment object to the cone mosaic object
theMosaic.os = theOuterSegment;


eyeMovementsPerOI = 0;  % no eye movements
[isomerizationRateSequence, photoCurrentSequence, responseTimeAxis] = ...
    coneMosaic_computeMultiOI(theConeMosaic, theOIsequence, eyeMovementsPerOI, 'currentFlag', true);
    

%% Plot the retinal photon rate at the center pixel and all wavelengths
figure(1); clf; hold on;
cmap = jet(size(referencePositionPhotons,1));
legends = {};
for waveIndex = 1:size(referencePositionPhotons,1)
    plot(stimulusTimeAxis,squeeze(referencePositionPhotons(waveIndex,:)), '-', 'Color', cmap(waveIndex,:), 'LineWidth', 1.5);
    legends{numel(legends)+1} = sprintf('%2.0f nm', oiWavelengthAxis(waveIndex));
end
legend(legends, 'Location', 'EastOutside');
xlabel('time (seconds');
ylabel('retinal photon rate at stimulus center');


%% Plot the isomerization rate sequence
figure(2); clf;
subplot(2,1,1);
plot(responseTimeAxis, squeeze(isomerizationRateSequence(1,1,:)), 'r-', 'LineWidth', 1.5);
hold on;
plot(responseTimeAxis, squeeze(isomerizationRateSequence(1,2,:)), 'g-', 'LineWidth', 1.5);
plot(responseTimeAxis, squeeze(isomerizationRateSequence(1,3,:)), 'b-', 'LineWidth', 1.5);
hold off;
ylabel('isomerization rate (R*/sec)');
xlabel('time (sec)');

subplot(2,1,2);
plot(responseTimeAxis, squeeze(photoCurrentSequence(1,1,:)), 'r-', 'LineWidth', 1.5);
hold on;
plot(responseTimeAxis, squeeze(photoCurrentSequence(1,2,:)), 'g-', 'LineWidth', 1.5);
plot(responseTimeAxis, squeeze(photoCurrentSequence(1,3,:)), 'b-', 'LineWidth', 1.5);
hold off;
ylabel('photocurrent');
xlabel('time (sec)');

return;




%% Set the stimulusSamplingInterval equal to the refresh cycle of a 120Hz display
stimulusSamplingInterval = 1/60;

%% Generate a cone mosaic object with cone
theConeMosaic = coneMosaic('size', [1 1]);
% Make it an L-cone (2)
theConeMosaic.pattern = 2;
% Make the integration time equal to the stimulus sampling interval
theConeMosaic.integrationTime = stimulusSamplingInterval;

%% Generate a linear outersegment object
theLinearOS = osLinear();
theLinearOS.noiseFlag = false;
% Set the outer-segment time step equal to the mosaic's integration time.
% This maybe ok for @osLinear, but definitely not for @osBiophys
% If we want a finer timerStep, we must correspondigly upsample the photonrate, before passing it to theMosaic.os.compute().
theLinearOS.timeStep = theConeMosaic.integrationTime;

%% Assign theLinearOS for to our mosaic
theMosaic.os = theLinearOS;

%% Generate a photon isomerization signal with a Gaussian-shaped modulation
% Time axis (-1.0 .. +1.0) seconds
timeAxis = -0.3:theConeMosaic.integrationTime:0.3;

meanPhotonRate = 3000;   % isomerizations / sec
peakModulationPhotonRate = 1000; % isomerizations / sec
isomerizationRateInPhotonsPerSec = meanPhotonRate + peakModulationPhotonRate * exp(-0.5*(timeAxis/(3*theConeMosaic.integrationTime)).^2);
isomerizationsPerTimeBin = isomerizationRateInPhotonsPerSec*theConeMosaic.integrationTime;
theMosaic.absorptions(1,1,:) = isomerizationRateInPhotonsPerSec;

% Compute totals
totalTime = numel(timeAxis)*theConeMosaic.integrationTime;
totalPhotonsFromRate = sum(isomerizationRateInPhotonsPerSec*theConeMosaic.integrationTime);
totalPhotons = sum(isomerizationsPerTimeBin);

% Compute the photocurrent
photocurrent = theMosaic.os.compute(theMosaic.absorptions(1,1,:),theConeMosaic.pattern);

hFig = figure(1); clf;
set(hFig, 'Position', [10 10 650 1280])

subplot(3,1,1);
plot(timeAxis, isomerizationsPerTimeBin, 'r*', 'LineWidth', 1.5);
hold on;
plot([0 0], [0 10000], 'k-');
set(gca, 'YLim', [0 100]);
set(gca, 'XLim', [timeAxis(1) timeAxis(end)]);
xlabel('time');
ylabel('photons / time bin');
title(sprintf('photon sequence, total photons: %g in %2.2f secs', totalPhotons, totalTime));

subplot(3,1,2);
stairs(timeAxis-theConeMosaic.integrationTime/2, isomerizationRateInPhotonsPerSec, 'r-', 'LineWidth', 1.5);
hold on
plot([0 0], [0 10000], 'k-');
set(gca, 'YLim', [0 5000]);
set(gca, 'XLim', [timeAxis(1) timeAxis(end)]);
xlabel('time');
ylabel('isomerization rate: photons / sec');
title(sprintf('isomerization rate, total photons: %g in %2.2f secs', totalPhotonsFromRate, totalTime));

subplot(3,1,3);
max(photocurrent(:))
stairs(timeAxis-theConeMosaic.integrationTime/2, squeeze(photocurrent(1,1,:)), 'r-', 'LineWidth', 1.5);
hold on
plot([0 0], [-85 0], 'k-');
set(gca, 'YLim', [-85 -60]);
set(gca, 'XLim', [timeAxis(1) timeAxis(end)]);
title('photocurrent');
drawnow;

end


% Supporting functions
function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
    uniformScene = sceneCreate('uniformd65');
    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);
    % set the radiance (in photons/steradian/m^2/nm)
    photonFlux = 1e16;
    uniformScene = sceneSet(uniformScene, 'photons', photonFlux*ones(64,64,numel(sceneGet(uniformScene, 'wave'))));
    % set desired luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end


% To be added to coneMosaic next to compute();
function [isomerizationRateSequence, photoCurrentSequence, responseTimeAxis] = coneMosaic_computeMultiOI(theConeMosaic, theOIsequence, eyeMovementsPerOI, varargin)

    class(theConeMosaic)
    % parse inputs
    p = inputParser;
    p.addRequired('theConeMosaic', @(x)isa(x, 'coneMosaic'));
    p.addRequired('theOIsequence', @iscell);
    p.addRequired('eyeMovementsPerOI', @isnumeric);
    p.addParameter('currentFlag', false, @islogical);
    p.parse(theConeMosaic, theOIsequence, eyeMovementsPerOI,varargin{:});
    currentFlag = p.Results.currentFlag;
    
    % Generate an expanded mosaic to deal with eye movements
    padRows = 64;
    padCols = 64;
    theExpandedMosaic = theConeMosaic.copy();
    theExpandedMosaic.pattern = zeros(theConeMosaic.rows+2*padRows, theConeMosaic.cols+2*padCols);
    theExpandedMosaic.noiseFlag = false;

    for oiIndex = 1:numel(theOIsequence)
        % Compute noise-free isomerizations at each cone location for the current optical image
        fullFrameIsomerizatios{oiIndex} = theExpandedMosaic.computeSingleFrame(theOIsequence{oiIndex},'FullLMS',true);
    end % oiIndex
    
    clear 'theExpandedMosaic'
    
    % Generate eye movements for the entire sequence of OIs.
    % If 'eyeMovementsPerOI' is <= 0, generate 1 eye movement/oi with zero amplitude movement
    if (eyeMovementsPerOI > 0)
        eyeMovementsNum = eyeMovementsPerOI*numel(theOIsequence);
        eyeMovementSequence = theConeMosaic.emGenSequence(eyeMovementsNum);
    else
        eyeMovementsPerOI = 1;
        eyeMovementsNum = eyeMovementsPerOI*numel(theOIsequence);
        eyeMovementSequence = 0*theConeMosaic.emGenSequence(eyeMovementsNum);
    end

    
    % Loop over the oi sequence
    for oiIndex = 1:numel(theOIsequence)
       % Apply current frame eye movements to the mosaic
       eyeMovementIndices = (oiIndex-1)*eyeMovementsPerOI + (1:eyeMovementsPerOI);
       theConeMosaic.emPositions = eyeMovementSequence(eyeMovementIndices,:);

        % Compute noise-free isomerizations for the current oi by applying eye movements during this oi
        theCurrentOpticalImageIsomerizations = ...
           theConeMosaic.applyEMPath(fullFrameIsomerizatios{oiIndex}, ...
            'padRows',padRows, 'padCols',padCols);

        % Add noise
        if (theConeMosaic.noiseFlag)
            theCurrentOpticalImageIsomerizations = theConeMosaic.photonNoise(theCurrentOpticalImageIsomerizations,'newNoise', true);
        end

        % Accumulate isomerizations by adding theCurrentOpticalImageIsomerizations 
        if (oiIndex == 1)
            isomerizationCountSequence = theCurrentOpticalImageIsomerizations;
        else
            isomerizationCountSequence = cat(3, isomerizationCountSequence, theCurrentOpticalImageIsomerizations);
        end
    end % for oiIndex
    
    % Compute isomerization rate
    totalTime = theConeMosaic.integrationTime * size(isomerizationCountSequence,3);
    responseTimeAxis = linspace(0,totalTime,size(isomerizationCountSequence,3));
    effectiveIntegrationTime = theConeMosaic.integrationTime/eyeMovementsPerOI;
    isomerizationRateSequence = isomerizationCountSequence / effectiveIntegrationTime;
    
    if (currentFlag)  
        if (abs(theConeMosaic.os.timeStep - effectiveIntegrationTime) > 1e-6)
            % Compute new effective integration time
            timeScale = effectiveIntegrationTime/theConeMosaic.os.timeStep;
            effectiveIntegrationTimeNew = effectiveIntegrationTime / timeScale;
            
            % Rescale the isomerizationRateSequence
            isomerizationRateSequence = isomerizationCountSequence * timeScale;
            resampleTimeAxis = 0:effectiveIntegrationTimeNew:totalTime;
            
            % Resample the isomerizationRateSequence
            fprintf('Resampling (factor: %2.2f) isomerizationRateSequence to match the time step of the OuterSegment\n', timeScale);
            isomerizationRateSequenceResampled = zeros(size(isomerizationRateSequence,1), size(isomerizationRateSequence,2), numel(resampleTimeAxis));
            for iRow = 1:size(isomerizationRateSequence,1)
                for iCol = 1:size(isomerizationRateSequence,2)
                    isomerizationRateSequenceResampled(iRow, iCol,:) = interp1(responseTimeAxis, squeeze(isomerizationRateSequence(iRow, iCol,:)), resampleTimeAxis, 'linear');
                end
            end
            
            fprintf('Computing photocurrent\n');
            photoCurrentSequenceResampled = theConeMosaic.os.osCompute(isomerizationRateSequenceResampled, theConeMosaic.pattern, 'append', false);
            clear 'isomerizationRateSequenceResampled'
            
            % Back to original sampling 
            fprintf('Resampling photocurrent to match the time base of the eye movements\n');
            photoCurrentSequence = zeros(size(isomerizationRateSequence,1), size(isomerizationRateSequence,2), numel(responseTimeAxis));
            for iRow = 1:size(isomerizationRateSequence,1)
                for iCol = 1:size(isomerizationRateSequence,2)
                    photoCurrentSequence(iRow, iCol,:) = interp1(resampleTimeAxis, squeeze(photoCurrentSequenceResampled(iRow, iCol,:)), responseTimeAxis, 'linear');
                end
            end
            clear 'photoCurrentSequenceResampled'
        else
            % match between os.timeStep and effectiveIntegrationTime
            photoCurrentSequence = theConeMosaic.os.osCompute(isomerizationRateSequence, theConeMosaic.pattern, 'append', false);
        end
        
    else
        photoCurrentSequence = [];
    end
end



% debug = false;
% if (debug)
% Get the optics
%     optics = oiGet(theOI, 'optics');
%     no OTF
%     optics = opticsSet(optics, 'model', 'diffraction limited');
%     set back the customized optics
%     theOI = oiSet(theOI,'optics', optics);
%     set the lens Transmittance to 1.0
%     theLens = oiGet(theOI,'lens');
%     lensTransmittanceOriginal = theLens.get('transmittance')
%     theLens = theLens.set('absorbance', zeros(size(theLens.wave)))
%     lensTransmittance = theLens.get('transmittance')
%     pause
%     
%     lensTransmittanceOriginal = theLens.transmittance;
%     zeroAbsorbance = zeros(size(lensTransmittanceOriginal));
%     theLens.absorbance =  zeroAbsorbance;
%     lensTransmittanceFull = theLens.transmittance
%     theOI = oiSet(theOI, 'lens', theLens);
% 
% 
% 
%     FOVs = [0.125 0.25 0.5 1 2 4 8 16 32 64];
%     meanLuminance = 300;
%     for k = 1:numel(FOVs)
%         uniformScene = uniformFieldSceneCreate(FOVs(k), meanLuminance);
%         meanSceneLuminance(k) = sceneGet(uniformScene, 'mean luminance');
% 
%         Compute the optical image
%         theOI = oiCompute(theOI, uniformScene);
%         meanRetinalIlluminance(k) = oiGet(theOI, 'mean illuminance');
%     end
% 
%     figure(1); clf;
%     subplot(1,2,1)
%     plot(FOVs, meanSceneLuminance, 'ks-');
%     xlabel('scene size (deg)');
%     title('scene luminance');
% 
%     subplot(1,2,2)
%     plot(FOVs, meanRetinalIlluminance, 'ks-');
%     xlabel('scene size (deg)');
%     title('retinal illuminance');
%     drawnow;
%     pause
% 
%     figure(101); clf;
%     subplot(1,2,1);
%     imshow(sceneGet(uniformScene, 'RGB'))
% 
% end %(debug)
