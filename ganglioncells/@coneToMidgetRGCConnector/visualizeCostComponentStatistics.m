function visualizeCostComponentStatistics(obj, axSpatial, axChromatic, theCostComponentsMatrix)

    totalCosts = theCostComponentsMatrix(:,1);
    spatialVarianceCosts = theCostComponentsMatrix(:,2);
    chromaticVarianceCosts = theCostComponentsMatrix(:,3);

    plotSpatialVarianceCostStatistics(axSpatial, spatialVarianceCosts);
    plotChromaticVarianceCostStatistics(axChromatic, chromaticVarianceCosts);

end

function plotSpatialVarianceCostStatistics(ax, spatialVarianceCost)
    medianSpatialVarianceCost = median((spatialVarianceCost(spatialVarianceCost>=0)));
    meanSpatialVarianceCost = mean((spatialVarianceCost(spatialVarianceCost>=0)));
    
    maxSpatialVariance = max(spatialVarianceCost);

    if (maxSpatialVariance < 0.5)
        spatialVarianceTicks = 0:0.05:0.5;
    elseif (maxSpatialVariance < 1.0)
        spatialVarianceTicks = 0:0.1:1.0;
    elseif (maxSpatialVariance < 2.5)
        spatialVarianceTicks = 0:0.25:2.5;
    elseif (maxSpatialVariance < 5)
        spatialVarianceTicks = 0:0.5:5;
    elseif (maxSpatialVariance < 10)
        spatialVarianceTicks = 0:1:10;
    elseif (maxSpatialVariance < 20)
        spatialVarianceTicks = 0:2:20;
    else
        spatialVarianceTicks = 0:5:maxSpatialVariance;
    end

    spatialVarianceTicks = 0:0.1:1;


    zeroInputRGCsNum = numel(find(spatialVarianceCost == -99));
    zeroInputRGCsPercentage = zeroInputRGCsNum/numel(spatialVarianceCost)*100;

    [countsSpatial,edgesSpatial] = histcounts(spatialVarianceCost,spatialVarianceTicks);
    countsSpatialPercentage = countsSpatial / numel(spatialVarianceCost)*100;

    maxCount = max(countsSpatialPercentage);
    maxCount = (floor(maxCount/10)+1)*10;
    maxCount = 100;

    width = 0.5*(edgesSpatial(2)-edgesSpatial(1));
    
    bar(ax,edgesSpatial(1:end-1), countsSpatialPercentage, 1, 'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0 0 0], 'LineWidth', 1.0);
    hold(ax, 'on');
    bar(ax, -0.05, zeroInputRGCsPercentage , width, 'FaceColor', [0 0 0]);
    plot(ax, meanSpatialVarianceCost*[1 1], [-0.05 maxCount], 'k--', 'LineWidth', 1.5);
    xlabel(ax,'spatial variance (\sigma) ( x RGC spacing)');
    ylabel(ax, 'percentage');
    xtickangle(ax,0)
    set(ax, 'XTick', spatialVarianceTicks, 'XLim', [-0.05 1.05]);
    set(ax, 'YLim', [0 maxCount], 'YTick', 0:10:100);
    set(ax, 'FontSize', 16)
    grid(ax,'on'); box(ax, 'off');
    title(ax, sprintf('median: %2.4f, mean: %2.4f', medianSpatialVarianceCost, meanSpatialVarianceCost));

end

function plotChromaticVarianceCostStatistics(ax, chromaticVarianceCost)
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

    bar(ax,1:(numel(edgesChroma)-1), countsChromaPercentage, 1, 'FaceColor',[1 0.75 0.75],'EdgeColor',[1 0 0], 'LineWidth', 1.0);
    hold(ax, 'on');
    bar(ax, -0.05, zeroInputRGCsPercentage, width, 'FaceColor', [0 0 0]);
    plot(ax, meanChromaticVarianceCostMappedToXScale*[1 1], [0 maxCount], 'k--', 'LineWidth', 1.5);
    set(ax,'XTick', 1:(numel(edgesChroma)-1), 'XLim', [-0.2 numel(edgesChroma)], 'XTickLabel', {'0', '1/5', '1/4', '1/3', '2/5', '1/2', '1'})
    set(ax, 'YLim', [0 maxCount], 'YTick', 0:10:100);
    xlabel(ax, 'chromatic variance');
    set(ax, 'FontSize', 16);
    grid(ax, 'on'); box(ax, 'off');
    title(ax, sprintf('mean: %2.4f', meanChromaticVarianceCost));

end
