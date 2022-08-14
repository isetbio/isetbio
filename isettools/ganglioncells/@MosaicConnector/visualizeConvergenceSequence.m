function visualizeConvergenceSequence(currentPass, costsMatrix, costsNames, ...
    netReassignments, maxPassesNum)

    % Visualize convergence
    hFig = figure(5050); clf;
    set(hFig, 'Position', [100 100 1000 800]);

    maxCostsDisplayed = min([3 size(costsMatrix,2)]);
    maxCost = max([0.01 min([1 max(costsMatrix(:))])]);

    for costIndex = 1:maxCostsDisplayed
        subplot(2,2,costIndex);
        theCostSequence = squeeze(costsMatrix(:,costIndex));
        plot(0:currentPass, theCostSequence, 'bo-', ...
            'MarkerFaceColor', [0.5 0.8 0.9], 'MarkerSize', 14, 'LineWidth', 1.5);
        set(gca, 'YLim', [0 maxCost], 'XLim', [0 currentPass+1], 'FontSize', 16);
        set(gca, 'YTick', 0:0.1:1.0, 'XTick', 0:2:maxPassesNum);
        grid on
        ylabel(sprintf('%s', costsNames{costIndex}));
        xlabel('iteration no');
    end

    subplot(2,2,4);
    bar(1:currentPass, netReassignments, 1, ...
        'FaceColor', [0.5 0.8 0.9], 'EdgeColor', [0 0 1], 'LineWidth', 1.0);
    set(gca, 'XLim', [0 currentPass+1], 'YLim', [0 max([1 max(netReassignments)])], ...
        'FontSize', 16);
    set(gca, 'XTick', 0:2:maxPassesNum);
    grid on
    xlabel('iteration no');
    ylabel('reassignments');
    drawnow;

end