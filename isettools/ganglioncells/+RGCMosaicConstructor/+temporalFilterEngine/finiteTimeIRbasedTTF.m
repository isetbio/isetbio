%
% RGCMosaicConstructor.temporalFilterEngine.finiteTimeIRbasedTTF
%

function [theTTFwithFiniteTimeIR, theFiniteTimeImpulseResponseData] = finiteTimeIRbasedTTF(...
        temporalFrequencySupportHz, ...
        theFullTimeTTF, ...
        finiteTimeIRdurationSeconds)

   
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


    % Generate time window that is zero < 0 and at t > finiteTimeIRdurationSeconds
    idxToZero = find(theImpulseResponseData.temporalSupportSeconds > finiteTimeIRdurationSeconds);
    idxToKeep = find(theImpulseResponseData.temporalSupportSeconds <= finiteTimeIRdurationSeconds);
    theWindow = theImpulseResponseData.amplitude * 0 + 1;
    theWindow(idxToZero) = 0;
    theTimeLimitedImpulseResponse = theImpulseResponseData.amplitude .* theWindow;

    % Compute the TTF of the time-limited IR
    theTTFwithFiniteTimeIR = RGCMosaicConstructor.temporalFilterEngine.oneSidedTTFfromTemporalImpulseResponse(theTimeLimitedImpulseResponse);
    

    % Return the finite time ImpulseResponseData
    theFiniteTimeImpulseResponseData = theImpulseResponseData;
    theFiniteTimeImpulseResponseData.temporalSupportSeconds = theImpulseResponseData.temporalSupportSeconds(idxToKeep);
    theFiniteTimeImpulseResponseData.amplitude = theImpulseResponseData.amplitude(idxToKeep);

    

    debug = false;
    if (debug)
        figure(55);clf
        plot(theImpulseResponseData.temporalSupportSeconds, theImpulseResponseData.amplitude*0, 'k-', 'LineWidth', 1.0);
        hold on
        plot(finiteTimeIRdurationSeconds*[1 1], [-1 1], 'k--');
        plot(theImpulseResponseData.temporalSupportSeconds, theImpulseResponseData.amplitude, 'r-', 'LineWidth', 3.0);
        plot(theImpulseResponseData.temporalSupportSeconds, theTimeLimitedImpulseResponse, 'b-', 'LineWidth', 1.5);
    
       
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


