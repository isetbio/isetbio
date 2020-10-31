function plotMovements(theAxesHandle, iteration, maxMovements, tolerance)
    if (~isempty(maxMovements))
        plot(theAxesHandle, 1:numel(maxMovements), maxMovements,'ks-'); 
        if (~isempty(tolerance))
            hold(theAxesHandle, 'on');
            plot(theAxesHandle, [1 iteration-1], tolerance*[1 1], 'r-');
        end
        set(theAxesHandle, 'YScale', 'log', 'YLim', [0.01 1], 'YTick', [0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1]);
        xlabel(theAxesHandle, 'iteration');
        ylabel(theAxesHandle, 'movement (microns)');
    end
end