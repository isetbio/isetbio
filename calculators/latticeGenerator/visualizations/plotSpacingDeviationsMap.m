function plotSpacingDeviationsMap(theAxesHandle, rfPositions, spacingDeviations, visualizedFOVMicrons)
    if (~isempty(spacingDeviations))
        % Sampling vector
        sampling = struct('minPos', 1, 'maxPos', 0.5*visualizedFOVMicrons, 'intervals', 100, 'scale', 'log');
        % Generate 2D map from scattered values
        [deviationMap, mapSupport] = mapFromScatteredPositions(rfPositions, spacingDeviations, sampling);

        contourf(theAxesHandle, mapSupport(:,:,1), mapSupport(:,:,2), deviationMap, 0:0.05:1.0, 'LineWidth', 0.1);
        ylabel(theAxesHandle, 'space (microns)');
        set(theAxesHandle, 'CLim', [0 1], 'ZLim', [0 1]);
        colormap(theAxesHandle,jet)
        title(theAxesHandle, sprintf('spacing deviations = %2.3f-%2.3f microns', min(deviationMap(:)), max(deviationMap(:))));
    end
    axis(theAxesHandle,'square')
end