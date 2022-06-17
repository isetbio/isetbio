function reportConeInputStatistics(...
            chromaticSpatialVarianceTradeoff, ...
            RGCRFinputs, RGCRFweights, ...
            theInputConeMosaic, ...
            localConeToRGCDensityRatio, varargin)
% Report statistics of the cone mosaic -> RGC mosaic
%
% Syntax:
%   RGCRFconnector.reportConeInputStatistics(...
%            chromaticSpatialVarianceTradeoff, ...
%            RGCRFinputs, RGCRFweights, ...
%            theInputConeMosaic, ...
%            localConeToRGCDensityRatio, varargin)
%
% Description:
%   Report statistics of the cone mosaic -> RGC mosaic
%
% Inputs:
%   chromaticSpatialVarianceTradeoff        
%   RGCRFinputs
%   RGCRFweights
%   theInputConeMosaic
%
% Outputs:
%    None
%
% Optional key/value pairs
%   'figureHandle'                  - The figure handle on which to render the figure
%   'axesHandles'                   - The axes handles on which to render the figure

%   
% History:
%   5/11/2022       NPC     Wrote it
%


    p = inputParser;
    p.addParameter('figureHandle', [], @(x)(isempty(x)||isa(x, 'handle')));
    p.addParameter('axesHandles', [], @(x)(isempty(x)||isa(x, 'cell')));
    p.parse(varargin{:});

    hFig = p.Results.figureHandle;
    ax = p.Results.axesHandles;

    % Unpack input
    allConePositions = theInputConeMosaic.coneRFpositionsMicrons;
    allConeSpacings = theInputConeMosaic.coneRFspacingsMicrons;
    allConeTypes = theInputConeMosaic.coneTypes;

    rgcsNum = numel(RGCRFinputs);

    totalCost = zeros(1, rgcsNum);
    spatialVarianceCost = totalCost;
    chromaticVarianceCost = totalCost;

    parfor iRGC = 1:rgcsNum 
        % Retrieve indices of input cones and their weights
        theConeInputs = RGCRFinputs{iRGC};
        theConeWeights = RGCRFweights{iRGC};

        [totalCost(iRGC), ~, spatialVarianceCost(iRGC), chromaticVarianceCost(iRGC)] = ...
            RGCRFconnector.costToMaintainInputCones(chromaticSpatialVarianceTradeoff, ...
            allConePositions(theConeInputs,:), ...
            allConeSpacings(theConeInputs), ...
            allConeTypes(theConeInputs), ...
            theConeWeights, localConeToRGCDensityRatio(iRGC));
    end

    meanSpatialVarianceCost = mean(spatialVarianceCost);
    meanChromaticVarianceCost = mean(chromaticVarianceCost);
    maxSpatialVariance = max(spatialVarianceCost);

    if (maxSpatialVariance < 5)
        spatialVarianceTicks = 0:0.5:5;
    elseif (maxSpatialVariance < 10)
        spatialVarianceTicks = 0:1:10;
    elseif (maxSpatialVariance < 20)
        spatialVarianceTicks = 0:2:20;
    else
        spatialVarianceTicks = 0:5:maxSpatialVariance;
    end


    [countsSpatial,edgesSpatial] = histcounts(spatialVarianceCost,spatialVarianceTicks);
    countsSpatial = countsSpatial / numel(allConeTypes)*100;

    chromaticVarianceTicks = [0 1/5 1/4 1/3 2/5 1/2 1];
    idx = find(chromaticVarianceTicks >= meanChromaticVarianceCost);
    upperIndex = idx(1)
    upperLimit = chromaticVarianceTicks(upperIndex);
    idx = find(chromaticVarianceTicks <= meanChromaticVarianceCost);
    lowerIndex = idx(end)
    lowerLimit = chromaticVarianceTicks(lowerIndex);
    m = (meanChromaticVarianceCost-lowerLimit)/(upperLimit-lowerLimit);
    meanChromaticVarianceCostMappedToXScale = lowerIndex + m *(upperIndex-lowerIndex);

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
    hold(axSpatial, 'on');
    plot(axSpatial, meanSpatialVarianceCost*[1 1], [0 maxCount], 'k--', 'LineWidth', 1.5);
    xlabel(axSpatial,'spatial variance, \sigma ( x cone spacing)');
    ylabel(axSpatial, 'percentage');
    xtickangle(axSpatial,0)
    set(axSpatial, 'XTick', spatialVarianceTicks, 'XLim', [-0.25 0.25+maxSpatialVariance]);
    set(axSpatial, 'YLim', [0 maxCount], 'YTick', 0:25:100);
    set(axSpatial, 'FontSize', 16)
    grid(axSpatial,'on'); box(axSpatial, 'off');
    title(axSpatial, sprintf('mean: %2.2f', meanSpatialVarianceCost));

    bar(axChromatic,1:(numel(edgesChroma)-1), countsChroma, 1, 'FaceColor',[1 0.75 0.75],'EdgeColor',[1 0 0], 'LineWidth', 1.0);
    hold(axChromatic, 'on');
    plot(axChromatic, meanChromaticVarianceCostMappedToXScale*[1 1], [0 maxCount], 'k--', 'LineWidth', 1.5);
    set(axChromatic,'XTick', 1:(numel(edgesChroma)-1), 'XLim', [0 numel(edgesChroma)], 'XTickLabel', {'0', '1/5', '1/4', '1/3', '2/5', '1/2', '1'})
    set(axChromatic, 'YLim', [0 maxCount], 'YTick', 0:25:100);
    xlabel(axChromatic, 'chromatic variance');
    set(axChromatic, 'FontSize', 16);
    grid(axChromatic, 'on'); box(axChromatic, 'off');
    title(axChromatic, sprintf('mean: %2.2f', meanChromaticVarianceCost));
end
