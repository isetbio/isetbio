function visualizeConvergenceSequence(currentPass, costsMatrix, costsNames, ...
    netReassignments, maxPassesNum, plotTitle, figNo)

    % Visualize convergence
    hFig = figure(figNo); clf;
    set(hFig, 'Position', [100 100 2000 800], 'Color', [1 1 1], 'Name', plotTitle);

    maxCostsDisplayed = 5;

    for costIndex = 1:maxCostsDisplayed
        theCostSequence = squeeze(costsMatrix(:,costIndex));
        minCost = min(theCostSequence(:));
        maxCost = max(theCostSequence(:));
        rangeOfCosts(costIndex) = maxCost-minCost;
    end

    maxCostRange = max([0.02 max(rangeOfCosts)]);
    dTick = max([0.01 ceil(maxCostRange/5*100)/100]);


    for costIndex = 1:maxCostsDisplayed
        
        theCostSequence = squeeze(costsMatrix(:,costIndex));

        minCost = min(theCostSequence(:));
        maxCost = max(theCostSequence(:));

        if (isnan(minCost))
            continue;
        end

        dCost = maxCostRange-(maxCost-minCost);

        minCost = minCost - dCost*0.5;
        maxCost = maxCost + dCost*0.5;

        if (minCost == maxCost)
            minCost = -1;
            maxCost = 1;
        end

        
        currentPass = numel(theCostSequence);
        if (currentPass < 10)
            currentPassTick = 1;
        elseif (currentPass < 20)
            currentPassTick = 2;
        elseif (currentPass < 50)
            currentPassTick = 5;
        elseif (currentPass < 100)
            currentPassTick = 10;
        elseif (currentPass < 200)
            currentPassTick = 20;
        elseif (currentPass < 500)
            currentPassTick = 50;
        else
            currentPassTick = 100;
        end

        subplot(2,3,costIndex);
        plot(1:currentPass, theCostSequence, 'bo-', ...
            'MarkerFaceColor', [0.5 0.8 0.9], 'MarkerSize', 14, 'LineWidth', 1.5);

        set(gca, 'YLim', [minCost-10*eps maxCost+10*eps], 'XLim', [-0.5 currentPass+0.5], 'FontSize', 20);
        set(gca, 'YTick', 0:dTick:maxCost, 'XTick', 0:currentPassTick:currentPass);
        grid on
        ylabel(sprintf('%s', costsNames{costIndex}));
        xlabel('iteration no');
    end

    subplot(2,3,6);
    currentPass = numel(netReassignments);
    bar(1:currentPass, netReassignments, 1, ...
        'FaceColor', [0.5 0.8 0.9], 'EdgeColor', [0 0 1], 'LineWidth', 1.0);
    set(gca, 'XLim', [-0.5 currentPass+0.5], 'YLim', [0 ceil(1.1*max([1 max(netReassignments)]))], ...
        'FontSize', 20);
    set(gca, 'XTick', 0:currentPassTick:currentPass);
    grid on
    xlabel('iteration no');
    ylabel(plotTitle);
    drawnow;

    
end