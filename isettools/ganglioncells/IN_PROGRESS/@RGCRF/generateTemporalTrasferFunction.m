function temporalTransferFunction = generateTemporalTrasferFunction(obj, p)
    % Equation (2) of Benardete & Kaplan '97:
    % "The receptive field of the primate P retinal ganglion cell, I: Linear dynamics"

    % This model is a cascade of 2 stages:
    % (1) A subtractive high-pass fiter and 
    % (2) A series of lowPassStagesNum number of low-pass filters
    omega = 2 * pi * obj.tfSupport;
    temporalTransferFunction = abs(p.gainImpulsesPerSecondPerUnitContrast) * exp(-1j * omega * p.initialDelayMsec/1e3) .* ...
        (1 - (p.gSubtractive) ./ (1 + 1j*omega*(p.tauHighPassMsec/1e3)) ) .* ...
        (                   1 ./ (1 + 1j*omega*(p.tauLowPassMsec/1e3))  ) .^ (p.lowPassStagesNum);
    
end