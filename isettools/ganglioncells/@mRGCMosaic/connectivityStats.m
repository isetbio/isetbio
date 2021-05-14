function stats = connectivityStats(obj, figNo)

    % Reset stats
    maxConeInputsPerRFcenter = 100;
    stats.rgcsNumWithNonZeroInputs = 0;
    stats.coneInputsPerRGC = zeros(1,maxConeInputsPerRFcenter);
    stats.eccSupport = logspace(log10(0.1), log10(30), 32);
    stats.eccTable = zeros(numel(stats.eccSupport), maxConeInputsPerRFcenter);
    stats.mixedInputRGCs = 0;
    stats.orphanRGCsNum = 0;
    
    stats.rgcsNum = size(obj.coneConnectivityMatrix,2);
    stats.conesNum = size(obj.coneConnectivityMatrix,1);
    
    for RGCindex = 1:stats.rgcsNum    
        % Find cones which are connected to this RGC
        connectivityVector = full(squeeze(obj.coneConnectivityMatrix(:, RGCindex)));
        inputConeIDs = find(connectivityVector > 0.01);
        inputsNum = numel(inputConeIDs);
        
        % Update stats
        if (inputsNum == 0)
            stats.orphanRGCsNum = stats.orphanRGCsNum + 1;
            continue;
        end
        
        stats.rgcsNumWithNonZeroInputs = stats.rgcsNumWithNonZeroInputs + 1;
        stats.coneInputsPerRGC(inputsNum) = stats.coneInputsPerRGC(inputsNum) + 1;
        
        inputTypesNum = numel(unique(obj.inputConeMosaic.coneTypes(inputConeIDs)));
        if (inputTypesNum>1)
            stats.mixedInputRGCs = stats.mixedInputRGCs + 1;
        end
        
        ecc = sqrt(sum((obj.rgcRFpositionsDegs(RGCindex,:)).^2,2));
        [~,eccIndex] = min(abs(stats.eccSupport-ecc));
        stats.eccTable(eccIndex, inputsNum) = stats.eccTable(eccIndex, inputsNum) + 1;
    end % RGCindex
    
    cellsNumPerEcc = sum(stats.eccTable,2);
    cellsNumPerEcc(cellsNumPerEcc==0) = nan;
    
    if (~isempty(figNo))
        figure(figNo); clf;
        plot(stats.eccSupport, 100*stats.eccTable(:,1)./cellsNumPerEcc, 'ks-', 'LineWidth', 1.5);
        hold on;
        plot(stats.eccSupport, 100*stats.eccTable(:,2)./cellsNumPerEcc, 'rs-', 'LineWidth', 1.5);
        plot(stats.eccSupport, 100*stats.eccTable(:,3)./cellsNumPerEcc, 'bs-', 'LineWidth', 1.5);
        plot(stats.eccSupport, 100*stats.eccTable(:,4)./cellsNumPerEcc, 'ms-', 'LineWidth', 1.5);
        plot(stats.eccSupport, 100*sum(stats.eccTable(:,5:end),2)./cellsNumPerEcc, 'gs-', 'LineWidth', 1.5);
        xlabel('eccentricity (degs)');
        ylabel('% of cells');
        legend({'1 cone input', '2 cone inputs', '3 cone inputs', '4 cone inputs', '5+ cone inputs'});
        set(gca, 'XScale', 'log', 'YScale', 'log','XLim', [0.1 35], 'XTick', [0.1 0.3 1 3 10 30], 'YLim', [0.01 105], 'YTick', [0.01 0.03 0.1 0.3 1 3 10 30 100]);
        set(gca, 'FontSize', 16);
    end
    
end
