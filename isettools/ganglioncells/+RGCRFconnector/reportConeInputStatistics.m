function reportConeInputStatistics(...
            chromaticSpatialVarianceTradeoff, ...
            RGCRFinputs, RGCRFweights, ...
            allConeRFpositions, ...
            allConeRFspacings, ...
            allConeTypes)

    rgcsNum = numel(RGCRFinputs);

    population = zeros(1,100);
    coneInputsNum = zeros(1, rgcsNum);
    totalCost = coneInputsNum;
    spatialVarianceCost = coneInputsNum;
    chromaticVarianceCost = coneInputsNum;

    for iRGC = 1:rgcsNum 
        % Retrieve indices of input cones and their weights
        theConeInputs = RGCRFinputs{iRGC};
        theConeWeights = RGCRFweights{iRGC};

        coneInputsNum(iRGC) = numel(theConeWeights);
        population(coneInputsNum(iRGC)) = population(coneInputsNum(iRGC)) + 1;
        
        
        [totalCost(iRGC), ~, spatialVarianceCost(iRGC), chromaticVarianceCost(iRGC)] = ...
            RGCRFconnector.costToMaintainInputCones(chromaticSpatialVarianceTradeoff, ...
            allConeRFpositions(theConeInputs,:), ...
            allConeRFspacings(theConeInputs), ...
            allConeTypes(theConeInputs), ...
            theConeWeights);
    end

    hFig = figure(111); clf;
    set(hFig, 'Color', [1 1 1]);
    subplot(1,2,1);
    plot(coneInputsNum, spatialVarianceCost, 'ro', 'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 12);
    grid on
    set(gca, 'XTick', 1:20, 'XLim', [1 10], 'YLim', [0 2], 'YTick', 0:0.2:2, 'FontSize', 16);
    
    xlabel('number of cone inputs');
    ylabel('spatial variance cost')

    subplot(1,2,2); hold on
    cMap = brewermap(20,'set1');

    for inputsNum = 1:max(coneInputsNum)
        idx = find(coneInputsNum == inputsNum);
        c = chromaticVarianceCost(idx);
        [N,edges] = histcounts(c, 0:0.25:1);
        stem(edges(1:end-1),N/sum(N), 'ks-', 'MarkerSize', 12, 'MarkerFaceColor', cMap(inputsNum,:));
    end
end

function render2Dhistogram(ax, coneInputsNum, chromaticVarianceCost)

    coneInputTicks = 1:10;
    chromaticVarianceTicks = 0:0.2:1;

    [counts,Xedges,Yedges] = histcounts2(coneInputsNum,chromaticVarianceCost,coneInputTicks,chromaticVarianceTicks)
    histogram2(ax,'XBinEdges',Xedges,'YBinEdges',Yedges,'BinCounts',counts, ...
        'FaceColor','flat','EdgeColor',[0 0 0]);
    cmap = brewermap(1024, 'blues');
    colormap(cmap);
    set(gca, 'Color', cmap(1,:));
    xlabel(ax,'number of cone inputs');
    ylabel(ax,'chromatic variance cost');

end
