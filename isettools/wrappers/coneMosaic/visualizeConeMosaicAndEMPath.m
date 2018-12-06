function visualizeConeMosaicAndEMPath(theConeMosaic, fixEMOBJ)

    figure(); clf;
    theConeMosaic.visualizeGrid(...
        'axesHandle', subplot('Position', [0.1 0.15 0.8 0.8]), ...
        'labelConeTypes', false, ...
        'backgroundColor', [1 1 1], ...
        'displayVisualDegs', true);
    set(gca, 'FontSize', 12, 'YTickLabel', {});
    xlabel('space (degs)')
    ylabel('');
    
    % Get the emPath in meters so we can plot it on the same scale as the cone
    % mosaic.
    emPathsMeters = fixEMOBJ.emPosMicrons * 1e-6;
    hold on;
    trialNo = 1;
    plot(emPathsMeters(trialNo,:,1),emPathsMeters(trialNo,:,2), 'w-', 'LineWidth', 3.0);
    plot(emPathsMeters(trialNo,:,1),emPathsMeters(trialNo,:,2), 'r-', 'LineWidth', 1.5);
    set(gca, 'FontSize', 20);
    
    figure(); clf
    nTimePoints = size(emPathsMeters,2);
    timeAxis = theConeMosaic.integrationTime*(1:nTimePoints);
    emPathDegs = emPathsMeters*1e6/theConeMosaic.micronsPerDegree;
    plot(timeAxis, emPathDegs(trialNo,:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(timeAxis, emPathDegs(trialNo,:,2), 'b-', 'LineWidth', 1.5);
    legend({'x-pos', 'y-pos'});
    set(gca, 'FontSize', 20);
    xlabel('time (seconds)');
    ylabel('space (degrees)');
    axis 'square'
end
