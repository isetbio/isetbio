function visualizeSFTuningOfUnderlyingModel(ax, iRGC, synthParams, responseGain, spatialFrequenciesCPDHR, maxResponse)
    
    rc = synthParams.centerCharacteristicRadiiDegs(iRGC);
    rs = synthParams.surroundCharacteristicRadiiDegs(iRGC);
    kc = synthParams.centerPeakSensitivities(iRGC);
    ks = synthParams.surroundPeakSensitivities(iRGC);
    
    centerResponseTuning = kc*pi*rc^2*exp(-(pi*rc*spatialFrequenciesCPDHR).^2);
    surroundResponseTuning = ks*pi*rs^2*exp(-(pi*rs*spatialFrequenciesCPDHR).^2);
    comboResponse = centerResponseTuning - surroundResponseTuning;
    
    gain = responseGain / max(comboResponse);
    plot(ax, spatialFrequenciesCPDHR, gain*comboResponse, '-', 'Color', [0.5 0.5 0.5], 'LineStyle', '-');
    plot(ax, spatialFrequenciesCPDHR, gain*centerResponseTuning, '-', 'Color', [0.5 0 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
    plot(ax, spatialFrequenciesCPDHR, gain*surroundResponseTuning, '-', 'Color', [0 0 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
end

