function reportConeInputStatistics(...
            chromaticSpatialVarianceTradeoff, ...
            RGCRFinputs, RGCRFweights, ...
            theInputConeMosaic, varargin)

    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandles', [], @(x)(isempty(x)||isa(x, 'cell')));
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandles;

    % Unpack input
    allConeRFpositions = theInputConeMosaic.coneRFpositionsMicrons;
    allConeRFspacings = theInputConeMosaic.coneRFspacingsMicrons;
    allConeTypes = theInputConeMosaic.coneTypes;

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

    spatialVarianceTicks = 0:0.1:1.0;
    [countsSpatial,edgesSpatial] = histcounts(spatialVarianceCost,spatialVarianceTicks);
    countsSpatial = countsSpatial / numel(allConeTypes)*100;

    chromaticVarianceTicks = [0 1/5 1/4 1/3 2/5 1/2 1];
    [countsChroma,edgesChroma] = histcounts(chromaticVarianceCost,chromaticVarianceTicks);
    countsChroma = countsChroma / numel(allConeTypes)*100;

    
    % Initialize figure
    if (isempty(ax))
        if (isempty(hFig))
            hFig = figure(); clf;
            set(hFig, 'Color', [1 1 1], 'Position', [10 10 900 450]);
        end
        axSpatial = subplot('Position', [0.07 0.15 0.42 0.8]);
        axChromatic = subplot('Position', [0.58 0.15 0.42 0.8]);
    else
        axSpatial = ax{1};
        axChromatic = ax{2};
    end

    maxCount = max([max(countsSpatial) max(countsChroma)]);
    maxCount = (floor(maxCount/10)+1)*10;
    % Spatial variance histogram
    bar(axSpatial,edgesSpatial(1:end-1), countsSpatial, 1, 'FaceColor',[0.85 0.85 0.85],'EdgeColor',[0 0 0], 'LineWidth', 1.0);
    xlabel(axSpatial,'spatial variance, \sigma ( x cone spacing)');
    ylabel(axSpatial, 'percentage');
    xtickangle(axSpatial,0)
    set(axSpatial, 'XTick', spatialVarianceTicks, 'XLim', [-0.1 1]);
    set(axSpatial, 'YLim', [0 maxCount], 'YTick', 0:25:100);
    set(axSpatial, 'FontSize', 16)
    grid(axSpatial,'on'); box(axSpatial, 'off');


    bar(axChromatic,1:(numel(edgesChroma)-1), countsChroma, 1, 'FaceColor',[1 0.75 0.75],'EdgeColor',[1 0 0], 'LineWidth', 1.0);
    set(axChromatic,'XTick', 1:(numel(edgesChroma)-1), 'XLim', [0 numel(edgesChroma)], 'XTickLabel', {'0', '1/5', '1/4', '1/3', '2/5', '1/2', '1'})
    set(axChromatic, 'YLim', [0 maxCount], 'YTick', 0:25:100);
    xlabel(axChromatic, 'chromatic variance');
    set(axChromatic, 'FontSize', 16);
    grid(axChromatic, 'on'); box(axChromatic, 'off');

end
