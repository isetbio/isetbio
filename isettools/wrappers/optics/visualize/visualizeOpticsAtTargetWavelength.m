function visualizeOpticsAtTargetWavelength(theOI, targetWavelength, ...
    spatialSupportArcMin, spatialFrequencySupportCyclePerDeg, varargin)

    p = inputParser;
    p.addParameter('extraOTFData', []);

    % Parse input
    p.parse(varargin{:});
    extraOTFData = p.Results.extraOTFData;
    
    hFig = figure(); clf;
    set(hFig, 'Color', [1 1 1]);
    ax1 = subplot(1,2,1);
    visualizePSF(theOI, targetWavelength, spatialSupportArcMin, ...
        'axesHandle', ax1);

    ax2 = subplot(1,2,2);
    visualizeOTF(theOI, targetWavelength, spatialFrequencySupportCyclePerDeg, ...
        'axesHandle', ax2, 'extraData', extraOTFData);
end
