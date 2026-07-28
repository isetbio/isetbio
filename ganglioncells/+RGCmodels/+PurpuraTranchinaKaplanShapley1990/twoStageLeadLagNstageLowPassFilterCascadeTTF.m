%
% RGCmodels.PurpuraTranchinaKaplanShapley1990.twoStageLeadLagNstageLowPassFilterCascadeTTF
%

function theTTF = twoStageLeadLagNstageLowPassFilterCascadeTTF(params, temporalFrequencySupportHz)
    
	% Get params
    gain = params(1);                         % Responsitivity at 0 Hz
    conductionDelaySeconds = params(2);       % 

    timeConstant1Seconds = params(3);         % tau 1
    timeConstant2Seconds = params(4);         % tau 2

    lowPassTimeConstant1Seconds = params(5);  % tau 3
    lowPassTimeConstant2Seconds = params(6);  % tau 4

    n1LowPassStagesNum = params(7);           % n1
    n2LowPassStagesNum = params(8);           % n2

    % Circular frequency in radians
    omega = 2 * pi * temporalFrequencySupportHz;

    % Delay filter
    theDelayFilterTTF = exp(-1i * omega * conductionDelaySeconds);

    % 2-stage lead-lag filter
    % Only when timeConstant1Seconds > timeConstant2Seconds
    % is the lead-lag filter a high-pass filter
    % (i.e., gain and phase increase with TF)
    theLeadLagFilterTTF = ((1 + 1i * omega * timeConstant1Seconds) ./ (1 + 1i * omega * timeConstant2Seconds) ) .^ (2);

    % N1-stage low-pass filter with lowPassTimeConstant1Seconds time constant 
    theLowPassFilter1TTF = (1 + 1i * omega * lowPassTimeConstant1Seconds) .^ (-n1LowPassStagesNum);

    % N2-stage low-pass filter with lowPassTimeConstant2Seconds time constant 
    if (n2LowPassStagesNum>0)
        theLowPassFilter2TTF = (1 + 1i * omega * lowPassTimeConstant2Seconds) .^ (-n2LowPassStagesNum);
    else
        theLowPassFilter2TTF = theLowPassFilter1TTF*0+1;
    end

    theTTF = gain * theDelayFilterTTF .* theLeadLagFilterTTF .* theLowPassFilter1TTF .* theLowPassFilter2TTF;
end


