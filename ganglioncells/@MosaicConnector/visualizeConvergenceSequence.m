function visualizeConvergenceSequence(currentPass, costsMatrix, costsNames, ...
    netReassignments, maxPassesNum, plotTitle, figNo)

    % Visualize convergence
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [100 100 1000 800], 'Color', [1 1 1], 'Name', plotTitle);

    maxCostsDisplayed = min([3 size(costsMatrix,2)]);
    % minCost = max([0 min(costsMatrix(:))-0.02]);
    % maxCost = max([minCost min([1 max(costsMatrix(:))])])+0.02;
    % 
    % minCost = min(costsMatrix(:));
    % maxCost = max(costsMatrix(:));
    % dCost = maxCost - minCost;
    % minCost = minCost - dCost*0.1;
    % maxCost = maxCost + dCost*0.1;

    for costIndex = 1:maxCostsDisplayed
        theCostSequence = squeeze(costsMatrix(:,costIndex));
        minCost = min(theCostSequence(:));
        maxCost = max(theCostSequence(:));
        rangeOfCosts(costIndex) = maxCost-minCost;
    end

    maxCostRange = max([0.02 max(rangeOfCosts)]);
    dTick = max([0.01 round(maxCostRange/5*100)/100]);

    for costIndex = 1:maxCostsDisplayed
        subplot(2,2,costIndex);
        theCostSequence = squeeze(costsMatrix(:,costIndex));

        minCost = min(theCostSequence(:));
        maxCost = max(theCostSequence(:));
        dCost = maxCostRange-(maxCost-minCost);

        minCost = minCost - dCost*0.5;
        maxCost = maxCost + dCost*0.5;

        if (minCost == maxCost)
            minCost = -1;
            maxCost = 1;
        end

        plot(0:currentPass, theCostSequence, 'bo-', ...
            'MarkerFaceColor', [0.5 0.8 0.9], 'MarkerSize', 14, 'LineWidth', 1.5);
        set(gca, 'YLim', [minCost maxCost], 'XLim', [-0.5 currentPass+0.5], 'FontSize', 20);
        set(gca, 'YTick', 0:dTick:1.0, 'XTick', 0:1:maxPassesNum);
        grid on
        ylabel(sprintf('%s', costsNames{costIndex}));
        xlabel('iteration no');
    end

    subplot(2,2,4);
    bar(1:currentPass, netReassignments, 1, ...
        'FaceColor', [0.5 0.8 0.9], 'EdgeColor', [0 0 1], 'LineWidth', 1.0);
    set(gca, 'XLim', [-0.5 currentPass+0.5], 'YLim', [0 ceil(1.1*max([1 max(netReassignments)]))], ...
        'FontSize', 20);
    set(gca, 'XTick', 0:1:maxPassesNum);
    grid on
    xlabel('iteration no');
    ylabel(plotTitle);
    drawnow;

    
end