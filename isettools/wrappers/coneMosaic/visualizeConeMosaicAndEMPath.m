function visualizeConeMosaicAndEMPath(theConeMosaic, fixEMOBJ)

    figure(); clf;
    ax = subplot('Position', [0.1 0.05 0.35 0.95]);
    theConeMosaic.visualizeGrid(...
        'axesHandle', ax, ...
        'labelConeTypes', false, ...
        'backgroundColor', [1 1 1], ...
        'ticksInVisualDegs', true);
    set(ax, 'FontSize', 18);
    xlabel(ax,'\it space (degs)')
    ylabel(ax,'\it space (degs)');
    
    % Get the emPath in meters so we can plot it on the same scale as the cone
    % mosaic.
    emPathsMeters = fixEMOBJ.emPosMicrons * 1e-6;
    hold on;
    trialNo = 1;
    plot(emPathsMeters(trialNo,:,1),emPathsMeters(trialNo,:,2), 'c-', 'LineWidth', 3.0);
    plot(emPathsMeters(trialNo,:,1),emPathsMeters(trialNo,:,2), 'b-', 'LineWidth', 1.5);

    
    ax = subplot('Position', [0.60 0.05 0.35 0.95]);
    nTimePoints = size(emPathsMeters,2);
    timeAxis = theConeMosaic.integrationTime*(1:nTimePoints);
    emPathDegs = emPathsMeters*1e6/theConeMosaic.micronsPerDegree;
    plot(timeAxis, emPathDegs(trialNo,:,1), 'r-', 'LineWidth', 1.5); hold on;
    plot(timeAxis, emPathDegs(trialNo,:,2), 'b-', 'LineWidth', 1.5);
    legend({'x-pos', 'y-pos'});
    set(gca, 'FontSize', 18);
    xlabel('\it time (seconds)');
    ylabel('\it space (degrees)');
    set(gca, 'YLim', 0.5*max(theConeMosaic.fov)*[-1 1]);
    axis 'square'
    set(gca, 'XTick', [0:0.1:1], 'YTick', [-0.3:0.1:0.3]);
    box on; grid on;
    
end
