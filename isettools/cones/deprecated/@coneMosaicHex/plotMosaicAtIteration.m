function plotMosaicAtIteration(obj, ax, iteration, radius)

    range = radius*[-1 1];
    % Get cone locs at desired iteration
    if (iteration == 0)
        obj.coneLocsHexGrid = obj.initialLattice;
    else
        obj.coneLocsHexGrid = squeeze(obj.latticeAdjustmentSteps(iteration,:,:));
    end
    idx = find(sqrt(sum(obj.coneLocsHexGrid.^2,2)) <= radius);
    obj.renderHexMesh(ax, obj.coneLocsHexGrid(idx,1), obj.coneLocsHexGrid(idx,2), [0.5 0.5 0.5], 'none', 0.8, 0.2, '-');
    hold on;
    plot(ax,obj.coneLocsHexGrid(idx,1), obj.coneLocsHexGrid(idx,2), 'ko', 'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5]);
    hold off
    set(ax, 'XLim', range, 'YLim', range, 'XTick', [], 'YTick', [], 'FontSize', 20, 'XColor', 'none', 'YColor', 'none');
    axis 'square'
    title(ax,sprintf('iteration %d', iteration));
    drawnow;
end