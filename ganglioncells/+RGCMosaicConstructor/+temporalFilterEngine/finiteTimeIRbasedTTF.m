%
% RGCMosaicConstructor.temporalFilterEngine.finiteTimeIRbasedTTF
%

function [theTTFwithFiniteTimeIR, theFiniteTimeImpulseResponseData] = finiteTimeIRbasedTTF(...
        temporalFrequencySupportHz, ...
        theFullTimeTTF, ...
        finiteTimeIRdurationSeconds, ...
        leftFiniteDurationWindowDurationSeconds, ...
        rightFiniteDurationWindowDurationSeconds)
   
    % Derive the full IR
    theImpulseResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                    theFullTimeTTF, temporalFrequencySupportHz, ...
                    'causal', false, ...
                    'upsample', 1);


    M = numel(temporalFrequencySupportHz);
    N  = 2 * (M - 1);
    samplingFrequency = 2 * temporalFrequencySupportHz(end);
    temporalSupportSeconds = (0:N-1) / samplingFrequency;
    diffs = abs(temporalSupportSeconds(:)-theImpulseResponseData.temporalSupportSeconds(:));
    assert(all(diffs<10*eps*ones(size(diffs))), 'temporal supports differ. how can this be?');


    if (1==2)
        % OLD WAY
        idxToZero = find(theImpulseResponseData.temporalSupportSeconds > finiteTimeIRdurationSeconds);
        % Generate time window that is zero < 0 and at t > finiteTimeIRdurationSeconds
        theWindow = theImpulseResponseData.amplitude * 0 + 1;
        theWindow(idxToZero) = 0;
        theTimeLimitedImpulseResponse = theImpulseResponseData.amplitude .* theWindow;
    end

    theImpulseResponseData.amplitude = windowedImpulseResponseFunction(...
        theImpulseResponseData.amplitude, theImpulseResponseData.temporalSupportSeconds, ...
        finiteTimeIRdurationSeconds, leftFiniteDurationWindowDurationSeconds, rightFiniteDurationWindowDurationSeconds);


    % Compute the TTF of the time-limited IR
    theTTFwithFiniteTimeIR = RGCMosaicConstructor.temporalFilterEngine.oneSidedTTFfromTemporalImpulseResponse(theImpulseResponseData.amplitude);
    

    % Return the finite time ImpulseResponseData
    idxToKeep = find(theImpulseResponseData.temporalSupportSeconds <= finiteTimeIRdurationSeconds);
    theFiniteTimeImpulseResponseData = theImpulseResponseData;
    theFiniteTimeImpulseResponseData.temporalSupportSeconds = theImpulseResponseData.temporalSupportSeconds(idxToKeep);
    theFiniteTimeImpulseResponseData.amplitude = theImpulseResponseData.amplitude(idxToKeep);

    

    debug = ~true;
    if (debug)
        figure(55);clf
        plot(theImpulseResponseData.temporalSupportSeconds, theImpulseResponseData.amplitude*0, 'k-', 'LineWidth', 1.0);
        hold on
        plot(finiteTimeIRdurationSeconds*[1 1], [-1 1], 'k--');
        plot(theImpulseResponseData.temporalSupportSeconds, theImpulseResponseData.amplitude, 'r-', 'LineWidth', 3.0);
        %plot(theImpulseResponseData.temporalSupportSeconds, theTimeLimitedImpulseResponse, 'b-', 'LineWidth', 1.5);
    
       
        figure(56);clf
        subplot(1,2,1);
        plot(theFiniteTimeImpulseResponseData.temporalSupportSeconds, theFiniteTimeImpulseResponseData.amplitude, 'r-', 'LineWidth', 1.5);
    
        subplot(1,2,2);
        plot(temporalFrequencySupportHz, abs(theFullTimeTTF), 'k-', 'LineWidth', 2);
        hold on
        plot(temporalFrequencySupportHz, abs(theTTFwithFiniteTimeIR), 'r-',  'LineWidth', 1.4);
        legend({'fulltime', 'finite time'})
        pause
    end

end


function theImpulseResponseAmplitude = windowedImpulseResponseFunction(theImpulseResponseAmplitude, temporalSupportSeconds, ...
    finiteTimeIRdurationSeconds, leftWindowDurationSeconds, rightWindowDurationSeconds)

    sigmaT = 0.5*max([10*eps leftWindowDurationSeconds]);
    theLeftTimeIndices = find(temporalSupportSeconds<=leftWindowDurationSeconds);
    leftWindowDurationSeconds = temporalSupportSeconds(theLeftTimeIndices(end));
    theLeftWindow = exp(-abs(temporalSupportSeconds-leftWindowDurationSeconds)/sigmaT);

    sigmaT = 0.5*max([10*eps rightWindowDurationSeconds]);

    theRightTimeIndices = find(temporalSupportSeconds>=finiteTimeIRdurationSeconds-rightWindowDurationSeconds);
    rightWindowDurationSeconds = temporalSupportSeconds(theRightTimeIndices(1));
    theRightWindow = exp(-abs(temporalSupportSeconds-rightWindowDurationSeconds)/sigmaT);


    debug = ~true;
    if (debug)
        figure(22);
        plot(temporalSupportSeconds, theImpulseResponseAmplitude, 'k-', 'LineWidth', 4);
        hold on;
        plot(temporalSupportSeconds,theRightWindow, 'r-', 'LineWidth', 1.5);
        plot(temporalSupportSeconds,theLeftWindow, 'b-', 'LineWidth', 1.5);
    end



    
    theWindowedAmplitude = theImpulseResponseAmplitude;
    theWindowedAmplitude(theLeftTimeIndices) = theWindowedAmplitude(theLeftTimeIndices) .* theLeftWindow(theLeftTimeIndices);
    theWindowedAmplitude(theRightTimeIndices) = theWindowedAmplitude(theRightTimeIndices) .* theRightWindow(theRightTimeIndices);
    
    theImpulseResponseAmplitude = theWindowedAmplitude;
    
   if (debug)
        plot(temporalSupportSeconds,theImpulseResponseAmplitude, 'g-', 'LineWidth', 1.5);
        pause
    end

end
