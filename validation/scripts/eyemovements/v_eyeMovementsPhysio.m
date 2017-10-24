function varargout = v_eyeMovementsPhysio(varargin)

    %% Initialize
    %ieInit;

    %% Reproduce identical random numbers
    rng('default'); rng(3);

    % mosaicType: choose from {'rect' 'hexRegDefault'  'hexRegCustomLambda' 'hexEccBased'}
    mosaicsExamined = {'hexRegDefault'}; % {'rect' 'hexRegDefault'  'hexRegCustomLambda' 'hexEccBased'};
    
    % Tremor only
    emFlag = [1 0 0];
    
    % Drift only
    %emFlag = [0 1 0];
    
    % Micro-saccades only
    %emFlag = [0 0 1];
    
    % Integration time (also em sample time): 0.1 millisecond
    integrationTimeSeconds = 0.5/1000;
    
    % Duration of each trial
    trialLengthSeconds = 10.0;
    
    % How many trials to compute
    nTrials = 200; % 500;
    
    for iMosaic = 1:numel(mosaicsExamined)
        mosaicType = mosaicsExamined{iMosaic};
        cm = generateMosaic(mosaicType, integrationTimeSeconds);
        [theEMpathsArcMin, emTimeAxis] = generateEMpaths(cm, emFlag, trialLengthSeconds, nTrials);
        if (all(emFlag == [0 0 1]))
            % Microsaccades only
            emData{iMosaic} = analyzeMicrosaccadeDynamics(theEMpathsArcMin, emTimeAxis);
            plotAnalyzedMicrosaccadeEMdata(emData, iMosaic, mosaicsExamined);
        elseif (all(emFlag == [1 0 0]))
            % Tremors only
            detrendingPolynomialOrder = 0;  % 0 (zero-mean), 1 (linear detrending), 2 (2nd order poly detrend) etc ... 
            emData{iMosaic} = analyzeTremorDynamics(theEMpathsArcMin, emTimeAxis, detrendingPolynomialOrder);
            plotAnalyzedTremorEMdata(emData, iMosaic, mosaicsExamined);
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
    for iTrial = 1:size(theEMpathsArcMin,1)
        plot(emTimeAxis, squeeze(theEMpathsArcMin(iTrial, :,1)), 'r-');
        plot(emTimeAxis, squeeze(theEMpathsArcMin(iTrial, :,2)), 'b-');
    end
    set(gca, 'XLim', [emTimeAxis(1) emTimeAxis(end)], 'YLim', 30*[-1 1]);
    title(sprintf('%s', mosaicsExamined{iMosaic}));
    xlabel('time (seconds)');
    ylabel('position (arc min)');
end


function emData = analyzeTremorDynamics(theEMpathsArcMin, emTimeAxis, detrendingPolynomialOrder)
    emData.theEMpathsArcMin = theEMpathsArcMin;
    emData.emTimeAxis = emTimeAxis;
    
    
    nPoints = size(theEMpathsArcMin,2);
    nFFT = 2^(1+round(log(nPoints)/log(2)));
    maxFreq = 1/(2*emTimeAxis(2)-emTimeAxis(1));
    deltaFreq = maxFreq / (nFFT/2); 
    frequencyAxis = (0:1:(nFFT-1)) * deltaFreq;
    
    amplitudeAxisArcMin = (0:2:60)/60;
    for iTrial = 1:size(theEMpathsArcMin,1)
        xPos = squeeze(theEMpathsArcMin(iTrial,:,1));
        yPos = squeeze(theEMpathsArcMin(iTrial,:,2));
        xDeltas = diff(xPos);
        yDeltas = diff(yPos);
        diffSignal = sqrt(xDeltas.^2 + yDeltas.^2);
        [d.amplitudes, d.amplitudeHistogram, d.times, d.durations, d.intervals] = analyzeMicroSaccadeDiffSignal(emTimeAxis,diffSignal, amplitudeAxisArcMin);
        d.diffSignal = diffSignal;
        
        [p,s,mu] = polyfit(emTimeAxis,xPos,detrendingPolynomialOrder);
        d.trendXpos = polyval(p,emTimeAxis,[],mu);
        d.detrendXpos = xPos - d.trendXpos;
    
        [p,s,mu] = polyfit(emTimeAxis,yPos,detrendingPolynomialOrder);
        d.trendYpos = polyval(p,emTimeAxis,[],mu);
        d.detrendYpos = yPos - d.trendYpos ;
        
        % perform FT analysis on the detrended positions
        xPosToAnalyze = d.detrendXpos;
        yPosToAnalyze = d.detrendYpos;
        
        ft = abs(fft(xPosToAnalyze, nFFT));
        d.spectrum.amplitudeXpos = ft(1:nFFT/2);
        ft = abs(fft(yPosToAnalyze, nFFT));
        d.spectrum.amplitudeYpos = ft(1:nFFT/2);
        d.spectrum.frequency = frequencyAxis(1:nFFT/2);
        emData.stats{iTrial} = d;
    end % iTrial
end

function emData = analyzeMicrosaccadeDynamics(theEMpathsArcMin, emTimeAxis)
    
    emData.theEMpathsArcMin = theEMpathsArcMin;
    emData.emTimeAxis = emTimeAxis;
    
    amplitudeAxisArcMin = 0:2:60;
    for iTrial = 1:size(theEMpathsArcMin,1)
        xPos = squeeze(theEMpathsArcMin(iTrial,:,1));
        yPos = squeeze(theEMpathsArcMin(iTrial,:,2));
        xDeltas = diff(xPos);
        yDeltas = diff(yPos);
        diffSignal = sqrt(xDeltas.^2 + yDeltas.^2);
        [d.amplitudes, d.amplitudeHistogram, d.times, d.durations, d.intervals] = analyzeMicroSaccadeDiffSignal(emTimeAxis,diffSignal, amplitudeAxisArcMin);
        d.diffSignal = diffSignal;
        emData.stats{iTrial} = d;
    end
end

function [amplitudes, amplitudeHistogram, times, durations, intervals] = analyzeMicroSaccadeDiffSignal(timeAxis, diffSignal, amplitudeAxisArcMin)
    amplitudes = [];
    amplitudeHistogram.count = [];
    amplitudeHistogram.amplitudeAxis = [];
    times = [];
    durations = [];
    intervals = [];
    
    idx = find(diffSignal > 0);
    if (numel(idx) == 0)
        return;
    end
    
    dt = (timeAxis(2)-timeAxis(1))*1000;
    minDeltaTime = 2/1000;
    movementIndex = 1; k = 1;
    while (k <= numel(idx))
        tSample = idx(k);
        durations(movementIndex) = 0;
        amplitudes(movementIndex) = 0;
        times(movementIndex) = timeAxis(tSample);
        stillOnSameMovement = true;
        while (k <= numel(idx)) && stillOnSameMovement
            durations(movementIndex) = durations(movementIndex) + dt;
            % accumulate path during the saccade trace
            amplitudes(movementIndex) = amplitudes(movementIndex) + diffSignal(idx(k));
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
    amplitudeHistogram.count = histcounts(amplitudes, amplitudeAxisArcMin);
    amplitudeHistogram.amplitudeAxis = amplitudeAxisArcMin(1:end-1);
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
    
function [theEMpathsArcMin, emTimeAxis] = generateEMpaths(cm, emFlag, trialLengthSeconds, nTrials)
    % Generate this many eye movements per trial
    eyeMovementsPerTrial = trialLengthSeconds/cm.integrationTime;

    %% Build an eye movement  and cone mosaic structure
    em = emCreate;
    em = emSet(em, 'em flag', emFlag);
    
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

function plotAnalyzedTremorEMdata(emData, iMosaic, mosaicsExamined)

    figure(2); clf;
    d = emData{iMosaic}; 

    for iTrial = 1:size(d.theEMpathsArcMin,1)
        % Stats for current trial
        stats = d.stats{iTrial};
        if (iTrial == 1)
            count = stats.amplitudeHistogram.count;
            ftAxis = stats.spectrum.frequency;
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
     
    meanOverAllTrialsFreq = mean(count,1);
    meanOverAllTrialsFreq = meanOverAllTrialsFreq/ max(meanOverAllTrialsFreq);
    
    meanOverAllTrialsSpectrumXpos = mean(ftAmpXpos,1);
    meanOverAllTrialsSpectrumYpos = mean(ftAmpYpos,1);
    
    % Display em and spectrum for one trial
    displayedEMtrial = 1;
    stats = d.stats{displayedEMtrial};
    
    % Plot current trial emData
    subplot(3,3,[1 2 3]);
    hold on
    plot(d.emTimeAxis, stats.trendXpos, '-', 'Color', [1 0.5 0.6], 'LineWidth', 5);
    plot(d.emTimeAxis, stats.trendYpos, '-', 'Color', [0.6 0.5 1.0], 'LineWidth', 5);
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), 'ro-');
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'bo-');
    
    grid on
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', max(max(abs(squeeze(d.theEMpathsArcMin(displayedEMtrial,:,:)))))*[-1 1]);
    set(gca, 'FontSize', 16);
    xlabel('time (seconds)');
    ylabel('position (min arc)');
    legend({'detrend X-pos', 'detrend Y-pos', 'eye X-pos', 'eye Y-pos'});
    
    % Plot diff signal for current trial and the detected tremors
    subplot(3,3,[4 5 6]); 
    stem(stats.times, stats.amplitudes*60, 'mo', 'MarkerFaceColor', [0.8 0.5 0.8]);
    hold on;
    plot(d.emTimeAxis(1:end-1), stats.diffSignal*60, 'k.');
    grid on
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', [0 max(stats.amplitudes*60)]);
    set(gca, 'FontSize', 16);
    xlabel('time (seconds)');
    ylabel(sprintf('tremor amplitude\n(arc sec)'));
    
    % Plot spectrum for current trial 
    subplot(3,3,7); 
    plot(stats.spectrum.frequency, stats.spectrum.amplitudeXpos, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(stats.spectrum.frequency, stats.spectrum.amplitudeYpos, 'b-', 'LineWidth', 1.5);
    grid on
    set(gca, 'XScale', 'log', 'YScale', 'log');
    set(gca, 'XLim', [stats.spectrum.frequency(1) stats.spectrum.frequency(end)], 'YLim', [0.1 max([max(stats.spectrum.amplitudeXpos) max(stats.spectrum.amplitudeYpos)])]);
    set(gca, 'FontSize', 16);
    xlabel('frequency (Hz)');
    ylabel(sprintf('amplitude'));
    
    % Plot mean spectrum across all trials
    subplot(3,3,8); 
    plot(stats.spectrum.frequency, meanOverAllTrialsSpectrumXpos, 'r-', 'LineWidth', 1.5);
    hold on
    plot(stats.spectrum.frequency, meanOverAllTrialsSpectrumYpos, 'b-', 'LineWidth', 1.5);
    grid on
    set(gca, 'XScale', 'log', 'YScale', 'log');
    set(gca, 'XLim', [stats.spectrum.frequency(1) stats.spectrum.frequency(end)], 'YLim', [0.1 max([max(meanOverAllTrialsSpectrumXpos) max(meanOverAllTrialsSpectrumXpos)])]);
    set(gca, 'FontSize', 16);
    xlabel('frequency (Hz)');
    ylabel(sprintf('amplitude'));
    
    % Plot mean spectrum across all trials
    subplot(3,3,9); 
    bar(stats.amplitudeHistogram.amplitudeAxis*60, meanOverAllTrialsFreq, 1, 'FaceColor', [0.6 0.6 0.6]);
    set(gca, 'XLim', [0 45]);
    grid on
    xlabel('tremor amplitude (arc sec)');
    ylabel(sprintf('normalized\nfreq. of occurence'));
    set(gca, 'FontSize', 16);
    axis 'square'
end

function plotAnalyzedMicrosaccadeEMdata(emData, iMosaic, mosaicsExamined)
    figure(2); clf;
    d = emData{iMosaic}; 
    % Accumulate stats over all trials
    for iTrial = 1:size(d.theEMpathsArcMin,1)
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
    end % iTrial
    
    meanOverAllTrialsFreq = mean(count,1);
    meanOverAllTrialsFreq = meanOverAllTrialsFreq/ max(meanOverAllTrialsFreq);
    
    % Display em and diff signal for one trial
    displayedEMtrial = 1;
    % Stats for current trial
    stats = d.stats{displayedEMtrial};
    
    % Plot current trial emData
    subplot(3,3,[1 2 3]);
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,1)), 'ro-');
    hold on
    plot(d.emTimeAxis, squeeze(d.theEMpathsArcMin(displayedEMtrial, :,2)), 'bo-');
    grid on
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', max(max(abs(squeeze(d.theEMpathsArcMin(displayedEMtrial,:,:)))))*[-1 1]);
    set(gca, 'FontSize', 16);
    xlabel('time (seconds)');
    ylabel('position (min arc)');
    legend({'eye X-pos', 'eye Y-pos'});
    
    % Plot diff signal for current trial and the detected saccades
    subplot(3,3,[4 5 6]); 
    stem(stats.times, stats.amplitudes, 'mo', 'MarkerFaceColor', [0.8 0.5 0.8]);
    hold on;
    plot(d.emTimeAxis(1:end-1), stats.diffSignal, 'k.');
    grid on
    set(gca, 'XLim', [d.emTimeAxis(1) d.emTimeAxis(end)], 'YLim', [0 max(stats.amplitudes)]);
    set(gca, 'FontSize', 16);
    xlabel('time (seconds)');
    ylabel(sprintf('micro-saccade amplitude\n(min arc)'));
    
   % Plot the amplitude histogram (mean over all trials)
    subplot(3,3,7);
    bar(stats.amplitudeHistogram.amplitudeAxis, meanOverAllTrialsFreq, 1, 'FaceColor', [0.6 0.6 0.6]);
    set(gca, 'XLim', [0 45]);
    grid on
    xlabel('microsaccade amplitude (arcmin)');
    ylabel(sprintf('normalized\nfreq. of occurence'));
    set(gca, 'FontSize', 16);
    axis 'square'
        
    % Plot scattegram of durations vs amplitudes
    subplot(3,3,8);
    plot(durations, amplitudes, 'ko');
    xlabel('duration (ms)');
    ylabel('amplitude (arc min)');
    set(gca, 'FontSize', 16);
    grid on
    axis 'square'
        
    % Plot scattegram of durations vs intervals
    subplot(3,3,9);    
    plot(durations, intervals, 'ko');
    xlabel('duration (ms)');
    ylabel('interval (ms)');
    set(gca, 'FontSize', 16);
    grid on
    axis 'square'
end
