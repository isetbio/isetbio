function visualizeSFTuningOfUnderlyingModel(ax, iRGC, synthParams, spatialFrequenciesCPDHR, maxResponse)
    rc = synthParams.centerRadiiDegs(iRGC);
    rs = synthParams.surroundRadiiDegs(iRGC);
    kc = synthParams.centerPeakSensitivities(iRGC);
    ks = synthParams.surroundPeakSensitivities(iRGC);
    
    centerResponseTuning = kc*pi*rc^2*exp(-(pi*rc*spatialFrequenciesCPDHR).^2);
    surroundResponseTuning = ks*pi*rs^2*exp(-(pi*rs*spatialFrequenciesCPDHR).^2);
    comboResponse = centerResponseTuning - surroundResponseTuning;
    
    g = maxResponse/max(comboResponse);
    comboResponse = g * comboResponse;
    centerResponseTuning = g * centerResponseTuning;
    surroundResponseTuning = g * surroundResponseTuning;
    
    plot(ax, spatialFrequenciesCPDHR, comboResponse, '-', 'Color', [0.5 0.5 0.5], 'LineStyle', '-');
    plot(ax, spatialFrequenciesCPDHR, centerResponseTuning, '-', 'Color', [0.5 0 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
    plot(ax, spatialFrequenciesCPDHR, surroundResponseTuning, '-', 'Color', [0 0 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
end

