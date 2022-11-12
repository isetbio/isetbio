function visualizeConvergenceSequence(currentPass, costsMatrix, costsNames, ...
    netReassignments, maxPassesNum, plotTitle, figNo)

    % Visualize convergence
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [100 100 1000 800], 'Name', plotTitle);

    maxCostsDisplayed = min([3 size(costsMatrix,2)]);
    minCost = max([0 min(costsMatrix(:))-0.02]);
    maxCost = max([minCost min([1 max(costsMatrix(:))])])+0.02;

    for costIndex = 1:maxCostsDisplayed
        subplot(2,2,costIndex);
        theCostSequence = squeeze(costsMatrix(:,costIndex));
        plot(0:currentPass, theCostSequence, 'bo-', ...
            'MarkerFaceColor', [0.5 0.8 0.9], 'MarkerSize', 14, 'LineWidth', 1.5);
        set(gca, 'YLim', [minCost maxCost], 'XLim', [0 currentPass+1], 'FontSize', 16);
        set(gca, 'YTick', 0:0.02:1.0, 'XTick', 0:2:maxPassesNum);
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
    ylabel(plotTitle);
    drawnow;

end