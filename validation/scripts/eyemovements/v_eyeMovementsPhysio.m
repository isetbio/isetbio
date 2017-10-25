function varargout = v_eyeMovementsPhysio(varargin)

    %% Initialize
    %ieInit;

    %% Reproduce identical random numbers
    rng('default'); rng(3);

    % mosaicType: choose from {'rect' 'hexRegDefault'  'hexRegCustomLambda' 'hexEccBased'}
    mosaicsExamined = {'hexRegDefault'}; % {'rect' 'hexRegDefault'  'hexRegCustomLambda' 'hexEccBased'};
    
    % eye movement type: choose from {'tremor', 'drift', 'microsaccade'}
    emType = 'microsaccade';
    
    % Integration time (also eye movement sample time)
    integrationTimeSeconds = 1/1000;
    
    % Duration of each trial
    trialLengthSeconds = 5.0;
    
    % How many trials to compute
    nTrials = 200;
    
    for iMosaic = 1:numel(mosaicsExamined)
        mosaicType = mosaicsExamined{iMosaic};
        cm = generateMosaic(mosaicType, integrationTimeSeconds);
        [theEMpathsArcMin, emTimeAxis] = generateEMpaths(cm, emType, trialLengthSeconds, nTrials);
        speedEstimationIntervalSeconds = 100/1000;
        switch emType
            case 'tremor'
                detrendingPolynomialOrder = 0;  % 0 (zero-mean), 1 (linear detrending), 2 (2nd order poly detrend) etc ... 
                emData{iMosaic} = analyzeTremorDynamics(theEMpathsArcMin, emTimeAxis, detrendingPolynomialOrder, speedEstimationIntervalSeconds);
                plotAnalyzedTremorEMdata(emData, iMosaic, mosaicsExamined);
            case 'drift'
                emData{iMosaic} = analyzeDriftDynamics(theEMpathsArcMin, emTimeAxis, speedEstimationIntervalSeconds);
                plotAnalyzedDriftEMdata(emData, iMosaic, mosaicsExamined);
            case 'microsaccade'
                emData{iMosaic} = analyzeMicrosaccadeDynamics(theEMpathsArcMin, emTimeAxis, speedEstimationIntervalSeconds);
                plotAnalyzedMicrosaccadeEMdata(emData, iMosaic, mosaicsExamined);
            otherwise
                error('Unknown em type: ''%s''.', emType);
        end
    end % iMosaic
    
    figure(1); clf;
    for iMosaic = 1:numel(mosaicsExamined)
        % Plot the raw eye movement data
        plotRawEMdata(emData, iMosaic, mosaicsExamined);
    end
end


function plotRawEMdata(emData, iMosaic, mosaicsExamined)
    d = emData{iMosaic};
    theEMpathsArcMin = d.theEMpathsArcMin;
    emTimeAxis = d.emTimeAxis;

    subplot(2,2,iMosaic);
    hold on;
    for iTrial = 1:1 % size(theEMpathsArcMin,1)
        plot(emTimeAxis, squeeze(theEMpathsArcMin(iTrial, :,1)), 'r-');
        plot(emTimeAxis, squeeze(theEMpathsArcMin(iTrial, :,2)), 'b-');
    end
    set(gca, 'XLim', [emTimeAxis(1) emTimeAxis(end)], 'YLim', 30*[-1 1]);
    title(sprintf('%s', mosaicsExamined{iMosaic}));
    xlabel('time (seconds)');
    ylabel('position (arc min)');
end


function emData = analyzeTremorDynamics(theEMpathsArcMin, emTimeAxis, detrendingPolynomialOrder, speedEstimationIntervalSeconds)
    emData.theEMpathsArcMin = theEMpathsArcMin;
    emData.emTimeAxis = emTimeAxis;
    
    
    nPoints = size(theEMpathsArcMin,2);
    nFFT = 2^(1+round(log(nPoints)/log(2)));
    deltaT = emTimeAxis(2)-emTimeAxis(1);
    maxFreq = 1/(2*deltaT);
    deltaFreq = maxFreq / (nFFT/2); 
    frequencyAxis = (0:1:(nFFT-1)) * deltaFreq;
    
    % tremor amplitude histogram bins
    amplitudeAxisArcMin = (0:2:60)/60;
    
    % analyze each trial separately
    for iTrial = 1:size(theEMpathsArcMin,1)
        xPos = squeeze(theEMpathsArcMin(iTrial,:,1));
        yPos = squeeze(theEMpathsArcMin(iTrial,:,2));
        [d.amplitudes, d.amplitudeHistogram, d.times, d.durations, d.intervals, d.diffSignal, d.decimatedSignal, d.speedSignal, d.speedHistogram] = ...
            analyzeDiffSignal(emTimeAxis, xPos, yPos, amplitudeAxisArcMin, speedEstimationIntervalSeconds);
        
        % detrend x/y position signals
        [detrendedXpos, d.detrendingXpos] = detrendSignals(emTimeAxis, xPos, detrendingPolynomialOrder);
        [detrendedYpos, d.detrendingYpos] = detrendSignals(emTimeAxis, yPos, detrendingPolynomialOrder);
        
        % perform FT analysis on the detrended positions        
        ft = abs(fft(detrendedXpos, nFFT));
        d.spectrum.amplitudeXpos = ft(1:nFFT/2) * 1/deltaT;
        ft = abs(fft(detrendedYpos, nFFT));
        d.spectrum.amplitudeYpos = ft(1:nFFT/2) * 1/deltaT;
        d.spectrum.frequency = frequencyAxis(1:nFFT/2);
        emData.stats{iTrial} = d;
    end % iTrial
end

function emData = analyzeDriftDynamics(theEMpathsArcMin, emTimeAxis, speedEstimationIntervalSeconds)

    emData.theEMpathsArcMin = theEMpathsArcMin;
    emData.emTimeAxis = emTimeAxis;
    
    % microsaccade amplitude histogram bins
    amplitudeAxisArcMin = 0:2:60;
    
    % analyze each trial separately
    for iTrial = 1:size(theEMpathsArcMin,1)
        xPos = squeeze(theEMpathsArcMin(iTrial,:,1));
        yPos = squeeze(theEMpathsArcMin(iTrial,:,2));
        [d.amplitudes, d.amplitudeHistogram, d.times, d.durations, d.intervals, d.diffSignal, d.decimatedSignal, d.speedSignal, d.speedHistogram] = ...
            analyzeDiffSignal(emTimeAxis, xPos, yPos, amplitudeAxisArcMin, speedEstimationIntervalSeconds);
        emData.stats{iTrial} = d;
    end % iTrial
end

function emData = analyzeMicrosaccadeDynamics(theEMpathsArcMin, emTimeAxis, speedEstimationIntervalSeconds)
    
    emData.theEMpathsArcMin = theEMpathsArcMin;
    emData.emTimeAxis = emTimeAxis;
    
    % microsaccade amplitude histogram bins
    amplitudeAxisArcMin = 0:2:60;
    
    % analyze each trial separately
    for iTrial = 1:size(theEMpathsArcMin,1)
        xPos = squeeze(theEMpathsArcMin(iTrial,:,1));
        yPos = squeeze(theEMpathsArcMin(iTrial,:,2));
        [d.amplitudes, d.amplitudeHistogram, d.times, d.durations, d.intervals, d.diffSignal, d.decimatedSignal, d.speedSignal, d.speedHistogram] = ...
            analyzeDiffSignal(emTimeAxis, xPos, yPos, amplitudeAxisArcMin, speedEstimationIntervalSeconds);
        emData.stats{iTrial} = d;
    end % iTrial
end

function [amplitudes, amplitudeHistogram, times, durations, intervals, diffSignal, decimatedSignal, speedSignal, speedHistogram] = analyzeDiffSignal(timeAxis, xPos, yPos, amplitudeAxisArcMin, speedEstimationIntervalSeconds)
    amplitudes = [];
    amplitudeHistogram.count = [];
    amplitudeHistogram.amplitudeAxis = [];
    times = [];
    durations = [];
    intervals = [];
    decimatedSignal = [];
    speedSignal = [];
    speedHistogram.count = [];
    speedHistogram.speedAxisMinutesArcPerSecond = [];
    
    % Compute diff signal
    xDeltas = diff(xPos);
    yDeltas = diff(yPos);
    diffSignal = sqrt(xDeltas.^2 + yDeltas.^2);
    idx = find(diffSignal > 0);
    if (numel(idx) == 0)
        return;
    end
        
    dt = (timeAxis(2)-timeAxis(1));

    % Compute end-points for each transition
    minDeltaTime = 2/1000;
    movementIndex = 1; k = 1;
    while (k <= numel(idx))
        tSample = idx(k);
        durations(movementIndex) = 0;
        amplitudes(movementIndex) = 0;
        times(movementIndex) = timeAxis(tSample);
        stillOnSameMovement = true;
        decimatedSignal(movementIndex,:) = [xPos(idx(k)) ; yPos(idx(k))];
        while (k <= numel(idx)) && stillOnSameMovement
            durations(movementIndex) = durations(movementIndex) + dt*1000;
            % accumulate path during the saccade trace
            amplitudes(movementIndex) = amplitudes(movementIndex) + diffSignal(idx(k));
            decimatedSignal(movementIndex,:) = [xPos(idx(k)+1) ; yPos(idx(k)+1)];
            if (k+1<=numel(idx))
                stillOnSameMovement = timeAxis(idx(k+1)) <= timeAxis(idx(k)+1) + minDeltaTime;  % (idx(k+1) == idx(k)+1);
            else
                stillOnSameMovement = false;
            end
            k = k + 1;
        end
        movementIndex = movementIndex + 1;
    end
    
    % Intervals in milliseconds
    intervals = [nan diff(times)*1000];
    
    % Amplitude histogram
    amplitudeHistogram.count = histcounts(amplitudes, amplitudeAxisArcMin);
    amplitudeHistogram.amplitudeAxis = amplitudeAxisArcMin(1:end-1) + 0.5*(amplitudeAxisArcMin(2)-amplitudeAxisArcMin(1));
    
    % Speed signal; compute speed in 100 millisecond intervals 
    halfBinWidth = round(speedEstimationIntervalSeconds/2/dt);
    actualspeedEstimationIntervalSeconds = (halfBinWidth*2+1)*dt
    
    speedXMinArcPerSecond = nan(1,numel(xPos)-2*halfBinWidth-1);
    speedYMinArcPerSecond = speedXMinArcPerSecond;
    for k = 1+halfBinWidth:numel(xPos)-halfBinWidth
        speedXMinArcPerSecond(k) = abs((xPos(k+halfBinWidth)-xPos(k-halfBinWidth)))/actualspeedEstimationIntervalSeconds;
        speedYMinArcPerSecond(k) = abs((yPos(k+halfBinWidth)-yPos(k-halfBinWidth)))/actualspeedEstimationIntervalSeconds;
    end
    speedMinArcPerSecond = sqrt(speedXMinArcPerSecond.^2 + speedYMinArcPerSecond.^2);
    speedSignal(:,1) = speedXMinArcPerSecond;
    speedSignal(:,2) = speedYMinArcPerSecond;
    
    % speed histogram
    speedAxisArcSecPerSecond = (0:10:1000);
    speedHistogram.count = histcounts(speedMinArcPerSecond*60, speedAxisArcSecPerSecond);
    speedHistogram.speedAxisArcSecPerSecond = speedAxisArcSecPerSecond(1:end-1) + 0.5*(speedAxisArcSecPerSecond(2)-speedAxisArcSecPerSecond(1));
end

function [detrendedSignal, detrendingSignal] = detrendSignals(timeAxis, signal, detrendingPolynomialOrder)    
    [p,s,mu] = polyfit(timeAxis,signal,detrendingPolynomialOrder);
    detrendingSignal = polyval(p,timeAxis,[],mu);
    detrendedSignal = signal - detrendingSignal;
end

function cm = generateMosaic(mosaicType, integrationTimeSeconds)
    % Common params
    % FOV in degrees
    fovDegs = 0.2;

    % Mosaic params
    switch mosaicType
        case 'rect'
            cm = coneMosaic();
            cm.setSizeToFOV(fovDegs);
            
        case 'hexRegDefault'
            resamplingFactor = 8;
            eccBasedConeDensity = false;
            cm = coneMosaicHex(resamplingFactor, ...
                'fovDegs', fovDegs, ...
                'eccBasedConeDensity', eccBasedConeDensity);
                
        case 'hexRegCustomLambda'
            resamplingFactor = 8;
            eccBasedConeDensity = false;
            customLamda = 3.0;
            cm = coneMosaicHex(resamplingFactor, ...
                'fovDegs', fovDegs, ...
                'eccBasedConeDensity', eccBasedConeDensity, ... 
                'customLambda', customLamda, ...
                'latticeAdjustmentPositionalToleranceF', 0.1, ...
                'latticeAdjustmentDelaunayToleranceF', 0.1);
        
        case 'hexEccBased'
            resamplingFactor = 8;
            eccBasedConeDensity = true;
            cm = coneMosaicHex(resamplingFactor, ...
                'fovDegs', fovDegs, ...
                'eccBasedConeDensity', eccBasedConeDensity, ... 
                'latticeAdjustmentPositionalToleranceF', 0.1, ...
                'latticeAdjustmentDelaunayToleranceF', 0.1);
        otherwise
            error('Incorrect mosaic type');
    end % switch
    
    % Set the integration time - also the time base for eye movements
    cm.integrationTime = integrationTimeSeconds;  
    
end 
    
function [theEMpathsArcMin, emTimeAxis] = generateEMpaths(cm, emType, trialLengthSeconds, nTrials)
    % Generate this many eye movements per trial
    eyeMovementsPerTrial = round(trialLengthSeconds/cm.integrationTime);

    %% Build an eye movement  and cone mosaic structure
    em = emCreate;
    
    switch emType
        case 'tremor'
            em = emSet(em, 'em flag', [1 0 0]);
        case 'drift'
            em = emSet(em, 'em flag', [0 1 0]);
        case 'microsaccade'
            em = emSet(em, 'em flag', [0 0 1]);
    end % switch
    
    %% Generate eye movement paths
    theEMpaths = zeros(nTrials, eyeMovementsPerTrial, 2);
    theEMpathsMicrons = theEMpaths;
    for iTrial= 1:nTrials
        [theEMpaths(iTrial, :,:), theEMpathsMicrons(iTrial, :,:)] = cm.emGenSequence(eyeMovementsPerTrial, 'em', em);
    end

    % Return em positions in arc min
    theEMpathsArcMin = theEMpathsMicrons / 300 * 60;
    emTimeAxis = cm.timeAxis;
end


function plotAnalyzedDriftEMdata(emData, iMosaic, mosaicsExamined)
    d = emData{iMosaic}; 
    nTrials = size(d.theEMpathsArcMin,1);
    
    for iTrial = 1:size(d.theEMpathsArcMin,1)
        % Stats for current trial
        stats = d.stats{iTrial};
        if (iTrial == 1)
            amplitudeCount = stats.amplitudeHistogram.count;
            speedCount = stats.speedHistogram.count;
        else
            if (~isempty(stats.amplitudeHistogram.count))
                amplitudeCount = cat(1, amplitudeCount, stats.amplitudeHistogram.count);
            end
            if (~isempty(stats.speedHistogram.count))
                speedCount = cat(1, speedCount, stats.speedHistogram.count);
            end
        end
    end % iTrial
    
    % Compute mean speed across trials
    meanOverAllTrialsSpeedFreq = mean(speedCount,1);
    meanOverAllTrialsSpeedFreq = meanOverAllTrialsSpeedFreq/ max(meanOverAllTrialsSpeedFreq);
    
    % Display em for one trial
    displayedEMtrial = 1;
    stats = d.stats{displayedEMtrial};
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1400 1200], 'Color', [1 1 1]);
    
    % Plot current trial emData
    subplot(3,3,[1 2]);
    hold on
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), 'r.-', 'MarkerSize', 14);
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'b.-', 'MarkerSize', 14);
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', max(max(abs(squeeze(d.theEMpathsArcMin(displayedEMtrial,:,:)))))*[-1 1]);
    legend({'eye X-pos', 'eye Y-pos'});
    addPlotLabels('time (seconds)', 'position (min arc)', 16, 'trial #1', 'grid on', '');
    
    subplot(3,3,3);
    plot(squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'k.-', 'MarkerSize', 14);
    XYrange = max(max(abs(d.theEMpathsArcMin(displayedEMtrial,:,:)))) * [-1 1];
    set(gca, 'XLim', XYrange, 'YLim', XYrange, 'XTick', -100:2:100, 'YTick', -100:2:100);
    addPlotLabels('position (min arc)', 'position (min arc)', 16, 'trial #1', 'grid on', 'axis square');
    
    subplot(3,3, [4 5]);
    plot(d.emTimeAxis(1:size(stats.speedSignal,1)), squeeze(stats.speedSignal(:,1))*60, 'r.-');
    hold on
    plot(d.emTimeAxis(1:size(stats.speedSignal,1)), squeeze(stats.speedSignal(:,2))*60, 'b.-');
    legend({'x-speed', 'y-speed'});
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', [0 sqrt(2.0)*max(stats.speedSignal(:)*60)]);
    addPlotLabels('time (seconds)', sprintf('speed \n(arc sec / second)'), 16, 'trial #1', 'grid on', '');
    
    subplot(3,3, [7 8]);
    plot(d.emTimeAxis(1:size(stats.speedSignal,1)), sqrt(squeeze(stats.speedSignal(:,1)).^2 + squeeze(stats.speedSignal(:,2)).^2)*60, 'k.-');
    legend({'speed'});
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', [0 sqrt(2.0)*max(stats.speedSignal(:)*60)]);
    addPlotLabels('time (seconds)', sprintf('speed \n(arc sec / second)'), 16, 'trial #1', 'grid on', '');
    
    % Plot speed histogram
    subplot(3,3,9); 
    bar(stats.speedHistogram.speedAxisArcSecPerSecond, meanOverAllTrialsSpeedFreq, 1, 'FaceColor', [0.6 0.6 0.6]);
    set(gca, 'XLim', [0 5*60], 'XTick', 0:60:1000);
    addPlotLabels('drift speed (arc sec / second)', sprintf('normalized\nfreq. of occurence'), 16, sprintf('trials: %d', nTrials), 'grid on', 'square axis');
    

end

function plotAnalyzedTremorEMdata(emData, iMosaic, mosaicsExamined)
 
    d = emData{iMosaic}; 
    nTrials = size(d.theEMpathsArcMin,1);
    
    for iTrial = 1:size(d.theEMpathsArcMin,1)
        % Stats for current trial
        stats = d.stats{iTrial};
        if (iTrial == 1)
            count = stats.amplitudeHistogram.count;
            ftAmpXpos = stats.spectrum.amplitudeXpos;
            ftAmpYpos = stats.spectrum.amplitudeYpos;
        else
            if (~isempty(stats.amplitudeHistogram.count))
                count = cat(1,count, stats.amplitudeHistogram.count);
            end
            ftAmpXpos = cat(1,ftAmpXpos, stats.spectrum.amplitudeXpos);
            ftAmpYpos = cat(1,ftAmpYpos, stats.spectrum.amplitudeYpos);
        end
    end % iTrial
     
    % Compute means across trials
    meanOverAllTrialsFreq = mean(count,1);
    meanOverAllTrialsFreq = meanOverAllTrialsFreq/ max(meanOverAllTrialsFreq);
    meanOverAllTrialsSpectrumXpos = mean(ftAmpXpos,1);
    meanOverAllTrialsSpectrumYpos = mean(ftAmpYpos,1);
    
    % Display em and spectrum for one trial
    displayedEMtrial = 1;
    stats = d.stats{displayedEMtrial};
    
    hFig = figure(2); clf;
    set(hFig, 'Position', [10 10 1400 1200], 'Color', [1 1 1]);
    
    % Plot current trial emData
    subplot(3,3,[1 2]);
    hold on
    plot(d.emTimeAxis, stats.detrendingXpos, '-', 'Color', [1 0.5 0.6], 'LineWidth', 5);
    plot(d.emTimeAxis, stats.detrendingYpos, '-', 'Color', [0.6 0.5 1.0], 'LineWidth', 5);
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), 'r.-', 'MarkerSize', 14);
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'b.-', 'MarkerSize', 14);
    legend({'detrending X-pos', 'detrending Y-pos', 'eye X-pos', 'eye Y-pos'});
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', max(max(abs(squeeze(d.theEMpathsArcMin(displayedEMtrial,:,:)))))*[-1 1]);
    addPlotLabels('time (seconds)', 'position (min arc)', 16, 'trial #1', 'grid on', '');
    
    subplot(3,3,3);
    plot(squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'k.-', 'MarkerSize', 14);
    XYrange = max(max(abs(d.theEMpathsArcMin(displayedEMtrial,:,:)))) * [-1 1];
    set(gca, 'XLim', XYrange, 'YLim', XYrange, 'XTick', -100:2:100, 'YTick', -100:2:100);
    addPlotLabels('position (min arc)', 'position (min arc)', 16, 'trial #1', 'grid on', 'axis square');
    
    % Plot diff signal for current trial and the detected tremors
    subplot(3,3,[4 5]); 
    stem(stats.times, stats.amplitudes*60, 'mo', 'MarkerFaceColor', [0.8 0.5 0.8], 'MarkerSize', 10);
    hold on;
    plot(d.emTimeAxis(1:end-1), stats.diffSignal*60, 'k.', 'MarkerSize', 14);
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', [0 max(stats.amplitudes*60)]);
    addPlotLabels('time (seconds)', sprintf('tremor amplitude\n(arc sec)'), 16, 'trial #1', 'grid on', '');
    
    % Plot spectrum for current trial 
    subplot(3,3,7); 
    plot(stats.spectrum.frequency, stats.spectrum.amplitudeXpos, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(stats.spectrum.frequency, stats.spectrum.amplitudeYpos, 'b-', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    set(gca, 'XLim', [stats.spectrum.frequency(1) stats.spectrum.frequency(end)], 'YLim', [0.1 max([max(stats.spectrum.amplitudeXpos) max(stats.spectrum.amplitudeYpos)])]);
    addPlotLabels('frequency (Hz)', 'amplitude', 16, 'trial #1', 'grid on', '');

    
    % Plot mean spectrum across all trials
    subplot(3,3,8); 
    plot(stats.spectrum.frequency, meanOverAllTrialsSpectrumXpos, 'r-', 'LineWidth', 1.5);
    hold on
    plot(stats.spectrum.frequency, meanOverAllTrialsSpectrumYpos, 'b-', 'LineWidth', 1.5);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    set(gca, 'XLim', [stats.spectrum.frequency(1) stats.spectrum.frequency(end)], 'YLim', [0.1 max([max(meanOverAllTrialsSpectrumXpos) max(meanOverAllTrialsSpectrumXpos)])]);
    addPlotLabels('frequency (Hz)', 'amplitude', 16, sprintf('trials: %d', nTrials), 'grid on', '');
    
    % Plot amplitude histogram
    subplot(3,3,9); 
    bar(stats.amplitudeHistogram.amplitudeAxis*60, meanOverAllTrialsFreq, 1, 'FaceColor', [0.6 0.6 0.6]);
    set(gca, 'XLim', [0 45]);
    addPlotLabels('tremor amplitude (arc sec)', sprintf('normalized\nfreq. of occurence'), 16, sprintf('trials: %d', nTrials), 'grid on', 'square axis');
end

function plotAnalyzedMicrosaccadeEMdata(emData, iMosaic, mosaicsExamined)
    figure(2); clf;
    d = emData{iMosaic}; 
    nTrials = size(d.theEMpathsArcMin,1);
    
    % Accumulate stats over all trials
    for iTrial = 1:nTrials
        % Stats for current trial
        stats = d.stats{iTrial};
        if (iTrial == 1)
            count = stats.amplitudeHistogram.count;
            durations = stats.durations;
            amplitudes = stats.amplitudes;
            intervals = stats.intervals;
        else
            if (~isempty(stats.amplitudeHistogram.count))
                count = cat(1,count, stats.amplitudeHistogram.count);
            end
            durations = cat(2, durations, stats.durations);
            amplitudes = cat(2, amplitudes, stats.amplitudes);
            intervals = cat(2, intervals , stats.intervals);
        end
        decimatedSignal{iTrial}.xyPos = stats.decimatedSignal;
    end % iTrial
    
    meanOverAllTrialsFreq = mean(count,1);
    meanOverAllTrialsFreq = meanOverAllTrialsFreq/ max(meanOverAllTrialsFreq);
    
    % Display em and diff signal for one trial
    displayedEMtrial = 1;
    % Stats for current trial
    stats = d.stats{displayedEMtrial};
    
    % Plot current trial emData
    subplot(3,3,[1 2]);
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), 'r.-', 'MarkerSize', 14);
    hold on
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'b.-', 'MarkerSize', 14);
    legend({'eye X-pos', 'eye Y-pos'});
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', max(max(abs(squeeze(d.theEMpathsArcMin(displayedEMtrial,:,:)))))*[-1 1]);
    addPlotLabels('time (seconds)', 'position (min arc)', 16, 'trial #1', 'grid on', '');
    
    subplot(3,3,3);
    plot(squeeze(decimatedSignal{displayedEMtrial}.xyPos(:,1)), squeeze(decimatedSignal{displayedEMtrial}.xyPos(:,2)),'mo-', 'MarkerFaceColor', [0.8 0.5 0.8], 'MarkerSize', 14);
    hold on;
    plot(squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'k.-', 'MarkerSize', 14);
    XYrange = max(max(abs(d.theEMpathsArcMin(displayedEMtrial,:,:)))) * [-1 1];
    set(gca, 'XLim', XYrange, 'YLim', XYrange, 'XTick', -100:2:100, 'YTick', -100:2:100);
    addPlotLabels('position (min arc)', 'position (min arc)', 16, 'trial #1', 'grid on', 'axis square');
    
    % Plot diff signal for current trial and the detected saccades
    subplot(3,3,[4 5]); 
    stem(stats.times, stats.amplitudes, 'mo', 'MarkerFaceColor', [0.8 0.5 0.8], 'MarkerSize', 14);
    hold on;
    plot(d.emTimeAxis(1:end-1), stats.diffSignal, 'k.', 'MarkerSize', 14);
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', [0 max(stats.amplitudes)]);
    addPlotLabels('time (seconds)', sprintf('micro-saccade amplitude\n(min arc)'), 16, 'trial #1', 'grid on', '');
    
    % Plot scattegram of durations vs amplitudes
    subplot(3,3,7);
    plot(durations, amplitudes, 'k.');
    set(gca, 'XLim', [0 15]);
    addPlotLabels('duration (ms)', 'amplitude (arc min)', 16, sprintf('%d trials', nTrials), 'grid on', 'square axis');
        
    % Plot scattegram of durations vs intervals
    subplot(3,3,8);    
    plot(durations, intervals, 'k.');
    set(gca, 'XLim', [0 15]);
    addPlotLabels('duration (ms)', 'interval (ms)', 16, sprintf('%d trials', nTrials), 'grid on', 'square axis');
    
    % Plot the amplitude histogram (mean over all trials)
    subplot(3,3,9);
    bar(stats.amplitudeHistogram.amplitudeAxis, meanOverAllTrialsFreq, 1, 'FaceColor', [0.6 0.6 0.6]);
    set(gca, 'XLim', [0 45]);
    addPlotLabels('microsaccade amplitude (arcmin)', sprintf('normalized\nfreq. of occurence'), 16, sprintf('%d trials', nTrials), 'grid on', 'square axis');
end

function  addPlotLabels(theXlabel, theYlabel, theFontSize, theTitle, theGrid, theAxis)
    xlabel(theXlabel, 'FontWeight', 'bold');
    ylabel(theYlabel, 'FontWeight', 'bold');
    set(gca, 'FontSize', theFontSize);
    if (strcmp(theGrid, 'grid on'))
        grid on
    else
        grid off;
    end
    if (strcmp(theAxis, 'square axis')) || (strcmp(theAxis, 'axis square'))
        axis 'square'
    end
    title(theTitle);
end
