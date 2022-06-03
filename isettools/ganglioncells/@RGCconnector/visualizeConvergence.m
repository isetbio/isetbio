function visualizeConvergence(currentPass, netTotalCostInitial, netTotalCost, ...
    netSpatialCostInitial, netSpatialCost, ...
    netChromaticCostInitial, netChromaticCost, ...
    netReassignments, maxPassesNum)

    % Visualize convergence
        hFig = figure(222); clf;
        set(hFig, 'Position', [100 100 1000 800]);
        subplot(2,2,1);
        plot(0:currentPass, [netTotalCostInitial netTotalCost], 'bo-', 'MarkerFaceColor', [0.5 0.8 0.9], 'MarkerSize', 14, 'LineWidth', 1.5);
        set(gca, 'YLim', [0 0.5], 'XLim', [0 currentPass+1], 'FontSize', 16);
        set(gca, 'YTick', 0:0.1:0.5, 'XTick', 0:2:maxPassesNum);
        grid on
        ylabel('net total cost');
        xlabel('iteration no');

        subplot(2,2,2);
        spatialCostSequence = [netSpatialCostInitial netSpatialCost];
        maxY = max([0.1 max(spatialCostSequence)]);
        plot(0:currentPass, spatialCostSequence, 'bo-', 'MarkerFaceColor', [0.5 0.8 0.9], 'MarkerSize', 14, 'LineWidth', 1.5);
        set(gca, 'YLim', [0 maxY], 'XLim', [0 currentPass+1], 'FontSize', 16);
        set(gca, 'YTick', 0:0.05:0.5, 'XTick', 0:2:maxPassesNum);
        grid on
        ylabel('net spatial cost');
        xlabel('iteration no');

        subplot(2,2,3);
        chromaCostSequence = [netChromaticCostInitial netChromaticCost];
        maxY = max([0.1 max(chromaCostSequence)]);
        plot(0:currentPass, chromaCostSequence, 'bo-', 'MarkerFaceColor', [0.5 0.8 0.9], 'MarkerSize', 14, 'LineWidth', 1.5);
        set(gca, 'YLim', [0 maxY], 'XLim', [0 currentPass+1], 'FontSize', 16);
        set(gca, 'YTick', 0:0.05:0.5, 'XTick', 0:2:maxPassesNum);
        grid on
        xlabel('iteration no');
        ylabel('net chromatic cost');

        subplot(2,2,4);
        plot(1:currentPass, netReassignments, 'bo-', 'MarkerSize', 14, 'MarkerFaceColor', [0.5 0.8 0.9], 'LineWidth', 1.5);
        set(gca, 'XLim', [0 currentPass+1], 'YLim',[0 max([1 max(netReassignments)])], 'FontSize', 16);
        set(gca, 'XTick', 0:2:maxPassesNum);
        grid on
        xlabel('iteration no');
        ylabel('cone reassignments');
        drawnow;

end