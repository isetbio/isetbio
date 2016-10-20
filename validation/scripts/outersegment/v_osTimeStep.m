function varargout = v_osTimeStep(varargin)
%
% Demonstrate simulations using three different timebases, one for stimuli (based on stimulus refresh rate), 
% one for absorptions and eye movements (based on coneMosaic.integrationTime), and a third one for 
% outer segment current computations (based on os.timeStep) 
% Also demonstrates usage of the computeForOISequence() method of @coneMosaic, which computes 
% absorptions and photocurrents for a sequence of sequentially presented optical images with eye movements.
%
% NPC, ISETBIO TEAM, 2016
%

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Init
ieInit;

% Reproduce identical random number
rng('default'); rng(1);


% Define number of response instances
instancesNum = 200;


% scene mean luminance
meanLuminance = 200;  
    
% Steady params
c0 = struct(...
    'mosaicSize', nan, ...                      % 1 L-, 1 M-, and 1 S-cone only
    'meanLuminance', meanLuminance, ...         % scene mean luminance
    'modulation', 1.0, ...                      % 100%  modulation against background
    'modulationRegion', 'FULL', ...             % modulate the central image (choose b/n 'FULL', and 'CENTER')
    'stimulusSamplingInterval',  0.12, ...      % 8.3- Hz stimulus refresh, e.g., 100 msec per optical image
    'osTimeStep', 1/1000, ...                   % 1 millisecond
    'integrationTime', nan, ...                 % we will vary this one
    'photonNoise', true, ...
    'osNoise', false);

    
% 60 ms integrationTime
stimulusConditionIndex = 1;
theCondition = c0;
theCondition.integrationTime = 60/1000;                  
c{stimulusConditionIndex} = theCondition;

[theConeMosaic60msec, theOIsequence60msec, ...
 oiTimeAxis60msec, absorptionsTimeAxis60msec, photoCurrentTimeAxis60msec, ...
 allInstancesAbsorptionsCountSequence60msec, ...
 allInstancesIsomerizationRateSequence60msec, ...
 allInstancesPhotoCurrents60msec] = runSimulation(c{stimulusConditionIndex}, instancesNum);  


% 30 ms integrationTime
stimulusConditionIndex = 2;
theCondition = c0;
theCondition.integrationTime = 30/1000;                  
c{stimulusConditionIndex} = theCondition;

[theConeMosaic30msec, theOIsequence30msec, ...
 oiTimeAxis30msec, absorptionsTimeAxis30msec, photoCurrentTimeAxis30msec, ...
 allInstancesAbsorptionsCountSequence30msec, ...
 allInstancesIsomerizationRateSequence30msec, ...
 allInstancesPhotoCurrents30msec] = runSimulation(c{stimulusConditionIndex}, instancesNum);  


% conditions validation data
UnitTest.validationData('condParams', c);

% 60 ms integration time validation data
UnitTest.validationData('oiTimeAxis60msec', oiTimeAxis60msec);
UnitTest.validationData('absorptionsTimeAxis60msec', absorptionsTimeAxis60msec);
UnitTest.validationData('photoCurrentTimeAxis60msec', photoCurrentTimeAxis60msec);
UnitTest.validationData('allInstancesAbsorptionsCountSequence60msec', allInstancesAbsorptionsCountSequence60msec);
UnitTest.validationData('allInstancesIsomerizationRateSequence60msec', allInstancesIsomerizationRateSequence60msec);
UnitTest.validationData('allInstancesPhotoCurrents60msec', allInstancesPhotoCurrents60msec, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'allInstancesPhotoCurrents60msec', 1e-7);

% 30 ms integration time validation data
UnitTest.validationData('oiTimeAxis30msec', oiTimeAxis30msec);
UnitTest.validationData('absorptionsTimeAxis30msec', absorptionsTimeAxis30msec);
UnitTest.validationData('photoCurrentTimeAxis30msec', photoCurrentTimeAxis30msec);
UnitTest.validationData('allInstancesAbsorptionsCountSequence30msec', allInstancesAbsorptionsCountSequence30msec);
UnitTest.validationData('allInstancesIsomerizationRateSequence30msec', allInstancesIsomerizationRateSequence30msec);
UnitTest.validationData('allInstancesPhotoCurrents30msec', allInstancesPhotoCurrents30msec, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'allInstancesPhotoCurrents30msec', 1e-7);

% Extra data: coneMosaics and oiSequences
UnitTest.extraData('theConeMosaic60msec', theConeMosaic60msec);
UnitTest.extraData('theOIsequence60msec', theOIsequence60msec);
UnitTest.extraData('theConeMosaic30msec', theConeMosaic30msec);
UnitTest.extraData('theOIsequence30msec', theOIsequence30msec);

% Plot results
if (runTimeParams.generatePlots)
    plotSNR(absorptionsTimeAxis60msec, oiTimeAxis60msec, photoCurrentTimeAxis60msec, ...
            allInstancesAbsorptionsCountSequence60msec, ...
            allInstancesPhotoCurrents60msec, ...
            1);

    plotSNR(absorptionsTimeAxis30msec, oiTimeAxis30msec, photoCurrentTimeAxis30msec, ...
            allInstancesAbsorptionsCountSequence30msec, ...
            allInstancesPhotoCurrents30msec, ...
            2);
end

end


function plotSNR(isomerizationsTimeAxis, oiTimeAxis, photocurrentTime, allInstancesIsomerizationsCount, allInstancesPhotoCurrents, figNo)
    
    % Compute isomerization means and stds
    isomerizationMeans = mean(allInstancesIsomerizationsCount, 4);
    isomerizationSTDs = std(allInstancesIsomerizationsCount, 0, 4);
    
    % Subtract first time point from all photocurrents
    timePoint = 1;
    photocurrentBaselineAtTimePoint1 = allInstancesPhotoCurrents(:,:,timePoint,:);
    allInstancesPhotoCurrents = bsxfun(@minus, allInstancesPhotoCurrents, reshape(photocurrentBaselineAtTimePoint1, [size(allInstancesPhotoCurrents,1) size(allInstancesPhotoCurrents,2) 1 size(allInstancesPhotoCurrents,4)])); 

    % compute photocurrent means and stds
    photocurrentMeansBaselineCorrected  = mean(allInstancesPhotoCurrents,4);
    photocurrentSTDsBaselineCorrected   = std(allInstancesPhotoCurrents, 0, 4);
    
    
    dt = isomerizationsTimeAxis(2)-isomerizationsTimeAxis(1);
    instancesNum = size(allInstancesIsomerizationsCount,4);
    % Plotting limits
    absorptionsFanoFactorLims = [0.0 10];
    photocurrentFanoFactorLims = [0.0 100];
    SNRlims = [0 60];
    photocurrentRange = [-5 50];  
    
    
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10+figNo*10 10 1800 1180], 'Color', [0 0 0]);
    
    colors = [1 0 0; 0 1.0 0; 0 0.8 1];
    coneNames = {'L-cone', 'M-cone', 'S-cone'};
    
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 3, ...
           'colsNum', 6, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.05, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.04);
    
    for coneType = 1:3
        mu = squeeze(isomerizationMeans(1,coneType,:));
        sigma = squeeze(isomerizationSTDs(1,coneType,:));
        % avoid divisions by zero
        sigma(mu == 0) = 1;
        
        variance = sigma.^2;
        isomerizationsInverseFanoFactor = mu ./ variance;
        isomerizationsSNR = mu ./sigma;
        
        mu = squeeze(photocurrentMeansBaselineCorrected(1,coneType,:));
        sigma = squeeze(photocurrentSTDsBaselineCorrected(1,coneType,:));
        % avoid divisions by near zero by making very small photocurrents = 0
        sigma(mu < 0.1) = 1;
        mu(mu < 0.1) = 0;
        
        variance = sigma.^2;
        photocurrentInverseFanoFactor = mu ./variance;
        photocurrentSNR = mu ./ sigma;
        
        maxIsomerizationCountForThisCone = max(max(max(squeeze(allInstancesIsomerizationsCount(:,coneType, :,:)))));
        minIsomerizationCountForThisCone = min(min(min(squeeze(allInstancesIsomerizationsCount(:,coneType, :,:)))));
        
        plotBackgroundColor = [0.1 0.1 0.1];
        
        % Absorption events
        subplot('Position', subplotPosVectors(coneType,1).v);
        hold on
        % Identify stimulus presentation times
        for k = 1:numel(oiTimeAxis)
            plot(oiTimeAxis(k)*[1 1], [minIsomerizationCountForThisCone maxIsomerizationCountForThisCone], 'k-', 'Color', [0.5 0.5 0.5]);
        end
        
        barOpacity = 0.1;
        for tIndex = 1:numel(isomerizationsTimeAxis)
            quantaAtThisTimeBin = squeeze(allInstancesIsomerizationsCount(:,coneType,tIndex,:));
            plot([isomerizationsTimeAxis(tIndex) isomerizationsTimeAxis(tIndex)+dt], [quantaAtThisTimeBin(:) quantaAtThisTimeBin(:)], '-', 'LineWidth', 1.5, 'Color', [colors(coneType,:) barOpacity]);
        end
        box on;
        set(gca, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Color', plotBackgroundColor, 'FontSize', 14, 'XLim', [isomerizationsTimeAxis(1) isomerizationsTimeAxis(end)+dt], 'YLim', [minIsomerizationCountForThisCone maxIsomerizationCountForThisCone]);
        if (coneType == 3)
            xlabel('time (seconds)', 'FontSize', 16);
        end
        ylabel(sprintf('%s isomerizations in %2.1f msec', coneNames{coneType}, 1000*dt), 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 1 0.3]);
        if (coneType == 1)
            title(sprintf('isomerizations\n(%d instances)', instancesNum), 'FontSize', 16, 'Color', [1 1 1]);
        end
        
        
        
        % Absorptions InverseFanoFactor
        subplot('Position', subplotPosVectors(coneType,2).v);
        hold on
        % Identify stimulus presentation times
        for k = 1:numel(oiTimeAxis)
            plot(oiTimeAxis(k)*[1 1], absorptionsFanoFactorLims, 'k-', 'Color', [0.5 0.5 0.5]);
        end
        % Identify the FanoFactor = 1 line
        plot([isomerizationsTimeAxis isomerizationsTimeAxis(end)+dt], [ones(size(isomerizationsTimeAxis)) 1], '--', 'Color', [0.8 0.8 0.3], 'LineWidth', 1.5);
        % Plot the time-varying Fano factor
        for tIndex = 1:numel(isomerizationsTimeAxis)
            if (tIndex < numel(isomerizationsTimeAxis))
                xx = [isomerizationsTimeAxis(tIndex) isomerizationsTimeAxis(tIndex)+dt isomerizationsTimeAxis(tIndex+1)];
                yy = [isomerizationsInverseFanoFactor(tIndex) * [1 1] isomerizationsInverseFanoFactor(tIndex+1)];
            else
                xx = [isomerizationsTimeAxis(tIndex) isomerizationsTimeAxis(tIndex)+dt];
                yy = isomerizationsInverseFanoFactor(tIndex) * [1 1];
            end
            plot(xx,yy, '-', 'Color', colors(coneType,:), 'LineWidth', 1.5);
        end
        box on;
        yTicks = [0 1 2 4 8];
        set(gca, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Color', plotBackgroundColor, 'FontSize', 14, 'XLim', [isomerizationsTimeAxis(1) isomerizationsTimeAxis(end)+dt], 'YLim', absorptionsFanoFactorLims, 'YScale', 'Linear', 'YTick', yTicks, 'YTickLabel', sprintf('%2.2f\n', yTicks));
        if (coneType == 3)
            xlabel('time (seconds)', 'FontSize', 16);
        end
        ylabel('quanta inverse Fano factor ( \mu/{\sigma}^2 )', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 1 0.3]);
        if (coneType == 1)
            title('Inv. Fano Factor ( \mu/{\sigma}^2 )', 'FontSize', 16, 'Color', [1 1 1]);
        end
        
        % Absorptions SNR
        subplot('Position', subplotPosVectors(coneType,3).v);
        hold on
        % Identify stimulus presentation times
        for k = 1:numel(oiTimeAxis)
            plot(oiTimeAxis(k)*[1 1], SNRlims, 'k-', 'Color', [0.5 0.5 0.5]);
        end
        % Plot the time-varying SNR
        for tIndex = 1:numel(isomerizationsTimeAxis)
            if (tIndex < numel(isomerizationsTimeAxis))
                xx = [isomerizationsTimeAxis(tIndex) isomerizationsTimeAxis(tIndex)+dt isomerizationsTimeAxis(tIndex+1)];
                yy = [isomerizationsSNR(tIndex) * [1 1] isomerizationsSNR(tIndex+1)];
            else
                xx = [isomerizationsTimeAxis(tIndex) isomerizationsTimeAxis(tIndex)+dt];
                yy = isomerizationsSNR(tIndex) * [1 1];
            end
            plot(xx,yy, '-', 'Color', colors(coneType,:), 'LineWidth', 1.5);
        end
        box on;
        set(gca, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Color', plotBackgroundColor, 'FontSize', 14, 'XLim', [isomerizationsTimeAxis(1) isomerizationsTimeAxis(end)+dt],  'YLim', SNRlims, 'YScale', 'Linear');
        if (coneType == 3)
            xlabel('time (seconds)', 'FontSize', 16);
        end
        ylabel('quanta SNR ( \mu/\sigma )', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 1 0.3]);
        if (coneType == 1)
            title('SNR ( \mu/\sigma )', 'FontSize', 16, 'Color', [1 1 1]);
        end
        
        
        % photocurrents
        subplot('Position', subplotPosVectors(coneType,4).v);
        % Identify stimulus presentation times
        hold on
        for k = 1:numel(oiTimeAxis)
            plot(oiTimeAxis(k)*[1 1], photocurrentRange, 'k-', 'Color', [0.5 0.5 0.5]);
        end
        
        % Plot photocurrents
        plot(photocurrentTime, squeeze(allInstancesPhotoCurrents(1, coneType, :, :)), 'LineWidth', 1.5, 'Color', [colors(coneType,:) barOpacity*2]);
        box on;
        set(gca, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Color', plotBackgroundColor, 'FontSize', 14, 'XLim', [isomerizationsTimeAxis(1) isomerizationsTimeAxis(end)+dt], 'YLim', photocurrentRange);
        if (coneType == 3)
            xlabel('time (seconds)', 'FontSize', 16);
        end
        ylabel('photocurrent (pA)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 1 0.3]);
        if (coneType == 1)
            title(sprintf('pcurrent, os.noiseFlag=false\n (%d instances)', instancesNum), 'FontSize', 16, 'Color', [1 1 1]);
        end

        % Photocurrent time-varying inverse FanoFactor
        
        subplot('Position', subplotPosVectors(coneType,5).v);
        hold on
        % Identify stimulus presentation times
        for k = 1:numel(oiTimeAxis)
            plot(oiTimeAxis(k)*[1 1], photocurrentFanoFactorLims, 'k-', 'Color', [0.5 0.5 0.5]);
        end
        % Identify the FanoFactor = 1 line
        plot([photocurrentTime photocurrentTime(end)], [ones(size(photocurrentTime)) 1], '--', 'Color', [0.8 0.8 0.3], 'LineWidth', 1.5);
        % Plot the time-varying Fano factor
        stairs(photocurrentTime, photocurrentInverseFanoFactor, 'Color', colors(coneType,:), 'LineWidth', 1.5);
        box on;
        yTicks = [1 5 10 30 50 100];
        set(gca, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Color', plotBackgroundColor, 'FontSize', 14, 'XLim', [isomerizationsTimeAxis(1) isomerizationsTimeAxis(end)+dt], 'YLim', photocurrentFanoFactorLims, 'YScale', 'Linear', 'YTick', yTicks, 'YTickLabel', sprintf('%2.2f\n', yTicks));
        if (coneType == 3)
            xlabel('time (seconds)', 'FontSize', 16);
        end
        ylabel('photocurrent Inverse Fano factor ( \mu/{\sigma}^2 )', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 1 0.3]);
        if (coneType == 1)
            title('Inv. Fano Factor ( \mu/{\sigma}^2 )', 'FontSize', 16, 'Color', [1 1 1]);
        end
        
        
        % photocurrents SNR
        subplot('Position', subplotPosVectors(coneType,6).v);
        hold on
        % Identify stimulus presentation times
        for k = 1:numel(oiTimeAxis)
            plot(oiTimeAxis(k)*[1 1], SNRlims, 'k-', 'Color', [0.5 0.5 0.5]);
        end
        % Plot the time-varying SNR
        stairs(photocurrentTime, photocurrentSNR, 'Color', colors(coneType,:), 'LineWidth', 1.5);
        box on;
        set(gca, 'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Color', plotBackgroundColor, 'FontSize', 14, 'XLim', [isomerizationsTimeAxis(1) isomerizationsTimeAxis(end)+dt],  'YLim', SNRlims, 'YScale', 'Linear');
        if (coneType == 3)
            xlabel('time (seconds)', 'FontSize', 16);
        end
        ylabel('photocurrent SNR ( \mu/\sigma )', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [1 1 0.3]);
        if (coneType == 1)
            title('SNR  ( \mu/\sigma )', 'FontSize', 16, 'Color', [1 1 1]);
        end
        
        
        drawnow;
    end
    
    %NicePlot.exportFigToPNG(sprintf('Fig%d.png', figNo), hFig, 300);
end

function [theConeMosaic, theOIsequence, ...
    oiTimeAxis, absorptionsTimeAxis, photoCurrentTimeAxis, ...
    allInstancesAbsorptionsCountSequence, ...
    allInstancesIsomerizationRateSequence, ...
    allInstancesPhotoCurrentSequence] = runSimulation(condData, instancesNum)

    mosaicSize = condData.mosaicSize;
    meanLuminance = condData.meanLuminance;
    modulation = condData.modulation;
    modulationRegion = condData.modulationRegion;
    stimulusSamplingInterval = condData.stimulusSamplingInterval;
    integrationTime = condData.integrationTime;
    osTimeStep = condData.osTimeStep;
    photonNoise = condData.photonNoise; 
    osNoise = condData.osNoise;
    
    % Define the time axis for the simulation (how much data we will generate)
    minTime = -0.84;
    maxTime = 0.6;
    oiTimeAxis = minTime:stimulusSamplingInterval:maxTime;
    stimulusRampTau = 0.165;

    % Generate a uniform field scene with desired mean luminance
    if (isnan(mosaicSize))
        FOV = 0.2;
    else
        FOV = max(mosaicSize);
    end
    theScene = uniformFieldSceneCreate(FOV, meanLuminance);

    % Generate optics
    noOptics = false;
    theOI = oiGenerate(noOptics);

    % Generate the sequence of optical images
    theOIsequence = oiSequenceGenerateForRampedSceneModulation(theScene, theOI, oiTimeAxis, stimulusRampTau, modulation, modulationRegion);

    % Generate the cone mosaic with eye movements for theOIsequence
    theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep, oiTimeAxis, numel(theOIsequence));

    % Make all movements 0.
    theConeMosaic.emPositions = theConeMosaic.emPositions * 0;
    
    % Compute instancesNum times (note: eyeMovements remain unchanged from instance to instance)
    for instanceIndex = 1:instancesNum
        [absorptionsCountSequence, absorptionsTimeAxis, photoCurrentSequence, photoCurrentTimeAxis] = ...
            theConeMosaic.computeForOISequence(theOIsequence, oiTimeAxis, ...
            'currentFlag', true, ...
            'newNoise', true ...
            );
        % Compute photon rate from photon count
        isomerizationRateSequence = absorptionsCountSequence / theConeMosaic.integrationTime;
    
        % Preallocate memory
        if (instanceIndex == 1)
            allInstancesAbsorptionsCountSequence = zeros([size(absorptionsCountSequence) instancesNum ]);
            allInstancesIsomerizationRateSequence = zeros([size(isomerizationRateSequence) instancesNum ]);
            allInstancesPhotoCurrentSequence = zeros([size(photoCurrentSequence) instancesNum ]);
        end
        
        allInstancesAbsorptionsCountSequence(:,:,:, instanceIndex) = absorptionsCountSequence;
        allInstancesIsomerizationRateSequence(:,:,:, instanceIndex) = isomerizationRateSequence;
        allInstancesPhotoCurrentSequence(:,:,:, instanceIndex) = photoCurrentSequence;
    end
end


% ------- Helper functions --------

function theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep, oiTimeAxis, opticalImageSequenceLength)
    % Default human mosaic
    theConeMosaic = coneMosaic;
    
    % Adjust size
    if isnan(mosaicSize)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicSize);
    end
    
    % Set the noise
    theConeMosaic.noiseFlag = photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = integrationTime;
    
    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = osNoise;
    
    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = osTimeStep;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;

    % Generate eye movement sequence for all oi's
    stimulusSamplingInterval = oiTimeAxis(2)-oiTimeAxis(1);
    eyeMovementsNumPerOpticalImage = stimulusSamplingInterval/theConeMosaic.integrationTime;
    eyeMovementsNum = round(eyeMovementsNumPerOpticalImage*opticalImageSequenceLength);
    
    if (eyeMovementsNum < 1)
        error('Less than 1 eye movement!!! \nStimulus sampling interval:%g Cone mosaic integration time: %g\n', stimulusSamplingInterval, theConeMosaic.integrationTime);
    else 
        fprintf('Optical image sequence contains %2.0f eye movements (%2.2f eye movements/oi)\n', eyeMovementsNum, eyeMovementsNumPerOpticalImage);
        theConeMosaic.emGenSequence(eyeMovementsNum);
    end
end


function theOIsequence = oiSequenceGenerateForRampedSceneModulation(theScene, theOI, oiTimeAxis, stimulusRampTau, modulation, modulationRegion)
    % Stimulus time ramp
    stimulusRamp = exp(-0.5*(oiTimeAxis/stimulusRampTau).^2);
    
    % Compute the optical image
    theOI = oiCompute(theOI, theScene);
    backgroundPhotons = oiGet(theOI, 'photons');
    
    % Preallocate cell array
    theOIsequence = cell(1,numel(oiTimeAxis));
    
    % Compute the seuence
    fprintf('Computing sequence of optical images\n');
    for stimFrameIndex = 1:numel(oiTimeAxis)
        if strcmp(modulationRegion, 'FULL')
            retinalPhotonsAtCurrentFrame = backgroundPhotons * (1.0 + modulation*stimulusRamp(stimFrameIndex));
        elseif strcmp(modulationRegion, 'CENTER')
            if (stimFrameIndex == 1)
                pos = oiGet(theOI, 'spatial support', 'microns');
                ecc = sqrt(squeeze(sum(pos.^2, 3)));
                idx = find(ecc < 0.5*max(pos(:)));
                [idx1, idx2] = ind2sub(size(ecc), idx);
            end
            retinalPhotonsAtCurrentFrame = backgroundPhotons;
            retinalPhotonsAtCurrentFrame(idx1, idx2, :) = retinalPhotonsAtCurrentFrame(idx1, idx2, :) * (1.0 + modulation*stimulusRamp(stimFrameIndex));
        else
            error('Unknown modulationRegion ''%s'', modulationRegion');
        end
        if (oiTimeAxis(stimFrameIndex) < -0.6)
            retinalPhotonsAtCurrentFrame = 0*retinalPhotonsAtCurrentFrame;
        end
        theOIsequence{stimFrameIndex} = oiSet(theOI, 'photons', retinalPhotonsAtCurrentFrame);
    end
end


function theOI = oiGenerate(noOptics)
    % Generate optics
    if (noOptics)
        theOI = oiCreate('diffraction limited');
        optics = oiGet(theOI,'optics');           
        optics = opticsSet(optics,'fnumber',0);
        optics = opticsSet(optics, 'off axis method', 'skip');
        theOI = oiSet(theOI,'optics', optics);
    else
        theOI = oiCreate('human');
    end
end


function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
    uniformScene = sceneCreate('uniformd65');
    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);
    % set some radiance (in photons/steradian/m^2/nm)
    photonFlux = 1e16;
    uniformScene = sceneSet(uniformScene, 'photons', photonFlux*ones(64,64,numel(sceneGet(uniformScene, 'wave'))));
    % adjust radiance according to desired  mean luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end


