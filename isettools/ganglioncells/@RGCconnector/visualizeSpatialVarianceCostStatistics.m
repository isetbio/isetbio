function visualizeSpatialVarianceCostStatistics(obj, axSpatial, spatialVarianceCost)
% Visualize spatial variance cost statistics

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

    zeroInputRGCsNum = numel(find(spatialVarianceCost == -99));
    zeroInputRGCsPercentage = zeroInputRGCsNum/numel(spatialVarianceCost)*100;

    [countsSpatial,edgesSpatial] = histcounts(spatialVarianceCost,spatialVarianceTicks);
    countsSpatialPercentage = countsSpatial / numel(spatialVarianceCost)*100;

    maxCount = max(countsSpatialPercentage);
    maxCount = (floor(maxCount/10)+1)*10;

    width = 0.5*(edgesSpatial(2)-edgesSpatial(1));
    
    bar(axSpatial,edgesSpatial(1:end-1), countsSpatialPercentage, 1, 'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0 0 0], 'LineWidth', 1.0);
    hold(axSpatial, 'on');
    bar(axSpatial, -0.05, zeroInputRGCsPercentage , width, 'FaceColor', [0 0 0]);
    plot(axSpatial, meanSpatialVarianceCost*[1 1], [-0.05 maxCount], 'k--', 'LineWidth', 1.5);
    xlabel(axSpatial,'spatial variance (\sigma) ( x RGC spacing)');
    ylabel(axSpatial, 'percentage');
    xtickangle(axSpatial,0)
    set(axSpatial, 'XTick', spatialVarianceTicks, 'XLim', [-0.05 0.25+maxSpatialVariance]);
    set(axSpatial, 'YLim', [0 maxCount], 'YTick', 0:10:100);
    set(axSpatial, 'FontSize', 16)
    grid(axSpatial,'on'); box(axSpatial, 'off');
    title(axSpatial, sprintf('median: %2.4f, mean: %2.4f', medianSpatialVarianceCost, meanSpatialVarianceCost));

end
