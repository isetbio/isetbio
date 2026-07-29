%
% RGCMosaicAnalyzer.visualize.dataAsBinnedHistogram(ax, theData, theBins, faceColor, edgeColor, faceAlpha, lineWidth, varargin)
%
function [theCounts, theBarPlotHandle] = dataAsBinnedHistogram(ax, theData, theBins, faceColor, edgeColor, faceAlpha, lineWidth, varargin)
	% Parse optional input
    p = inputParser;
    p.addParameter('relativeHeight', 1, @isscalar);
    p.addParameter('lineStyle', '-', @ischar);
    p.addParameter('binWidth', 1.0, @(x)(isempty(x)||isscalar(x)));
    p.addParameter('plotWithNegativePolarity', false, @islogical);

    p.parse(varargin{:});
    relativeHeight = p.Results.relativeHeight;
    lineStyle = p.Results.lineStyle;
    binWidth = p.Results.binWidth;
    plotWithNegativePolarity = p.Results.plotWithNegativePolarity;

	theCounts = [0 histcounts(theData(:), theBins)];
	theFrequency = theCounts/sum(theCounts(:))*relativeHeight;

    theBarPlotHandle = RGCMosaicAnalyzer.visualize.shadedHistogram(ax, theBins, theFrequency, binWidth, plotWithNegativePolarity, ...
        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);

end