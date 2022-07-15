function visualizeChromaticVarianceCostStatistics(obj, axChromatic, chromaticVarianceCost)
% Visualize chromatic variance cost

    meanChromaticVarianceCost = mean(chromaticVarianceCost(chromaticVarianceCost>=0));

    chromaticVarianceTicks = [0 1/5 1/4 1/3 2/5 1/2 1];
    idx = find(chromaticVarianceTicks >= meanChromaticVarianceCost);
    upperIndex = idx(1);
    upperLimit = chromaticVarianceTicks(upperIndex);
    idx = find(chromaticVarianceTicks <= meanChromaticVarianceCost);
    lowerIndex = idx(end);
    lowerLimit = chromaticVarianceTicks(lowerIndex);
    m = (meanChromaticVarianceCost-lowerLimit)/(upperLimit-lowerLimit);
    meanChromaticVarianceCostMappedToXScale = lowerIndex + m *(upperIndex-lowerIndex);

    zeroInputRGCsNum = numel(find(chromaticVarianceCost == -99));
    zeroInputRGCsPercentage = zeroInputRGCsNum/numel(chromaticVarianceCost)*100;

    [countsChroma,edgesChroma] = histcounts(chromaticVarianceCost,chromaticVarianceTicks);
    countsChromaPercentage = countsChroma / numel(chromaticVarianceCost)*100;

    maxCount = max(countsChromaPercentage);
    maxCount = (floor(maxCount/10)+1)*10;
    maxCount = 100;
    
    width = 0.2;

    bar(axChromatic,1:(numel(edgesChroma)-1), countsChromaPercentage, 1, 'FaceColor',[1 0.75 0.75],'EdgeColor',[1 0 0], 'LineWidth', 1.0);
    hold(axChromatic, 'on');
    bar(axChromatic, -0.05, zeroInputRGCsPercentage, width, 'FaceColor', [0 0 0]);
    plot(axChromatic, meanChromaticVarianceCostMappedToXScale*[1 1], [0 maxCount], 'k--', 'LineWidth', 1.5);
    set(axChromatic,'XTick', 1:(numel(edgesChroma)-1), 'XLim', [-0.2 numel(edgesChroma)], 'XTickLabel', {'0', '1/5', '1/4', '1/3', '2/5', '1/2', '1'})
    set(axChromatic, 'YLim', [0 maxCount], 'YTick', 0:10:100);
    xlabel(axChromatic, 'chromatic variance');
    set(axChromatic, 'FontSize', 16);
    grid(axChromatic, 'on'); box(axChromatic, 'off');
    title(axChromatic, sprintf('mean: %2.4f', meanChromaticVarianceCost));

end