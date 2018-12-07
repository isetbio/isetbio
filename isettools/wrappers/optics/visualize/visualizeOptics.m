function visualizeOptics(theOI, targetWavelength, ...
    spatialSupportArcMin, spatialFrequencySupportCyclePerDeg)

    figure(); clf;
    ax1 = subplot(1,2,1);
    visualizePSF(theOI, targetWavelength, spatialSupportArcMin, ...
        'axesHandle', ax1);

    ax2 = subplot(1,2,2);
    visualizeOTF(theOI, targetWavelength, spatialFrequencySupportCyclePerDeg, ...
        'axesHandle', ax2);
end
