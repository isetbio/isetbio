%
% RGCmodels.BenardeteKaplan1997.oneStageHighPassNstageLowPassFilterCascadeTTF
%

function theTTF = oneStageHighPassNstageLowPassFilterCascadeTTF(params, temporalFrequencySupportHz)

	% Get params
    gain = params(1);                         % Responsitivity at 0 Hz
    conductionDelaySeconds = params(2);       % D
    highPassGain = params(3);                 % Hs
    highPassTimeConstantSeconds = params(4);  % Tau_s (Benardete and Kaplan (1992a) varied this for different contrast levels)
    lowPassTimeConstantSeconds = params(5);   % Tau_l
    nLowPassStagesNum = params(6);            % Nl
    nHighPassStagesNum = params(7);           % always 1 in Benardete & Kaplan (1992a)

    % Circular frequency in radians
    omega = 2 * pi * temporalFrequencySupportHz;

    % Delay filter
    theDelayFilterTTF = exp(-1i * omega * conductionDelaySeconds);

    % 1-stage high-pass filter
    theHighPassFilterTTF = 1 - highPassGain * (1 + 1i * omega * highPassTimeConstantSeconds) .^ (-nHighPassStagesNum);

    % N-stge low-pass filter
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-nLowPassStagesNum);

    theTTF = gain * theDelayFilterTTF .* theHighPassFilterTTF .* theLowPassFilterTTF;
end

