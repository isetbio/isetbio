function theBarPlotHandle = shadedHistogram(ax, theBins, theFrequency, theBinWidth, plotWithNegativePolarity, ...
	faceColor, edgeColor, faceAlpha, lineWidth, lineStyle, varargin)
%
% RGCMosaicAnalyzer.visualize.shadedHistogram(ax,theBins, theFrequency, theBinWidth, plotWithNegativePolarity,...
%  faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
%

	p = inputParser;
	p.addParameter('flipXY', false, @islogical);
	p.parse(varargin{:});
	flipXY = p.Results.flipXY;

	[xx, yy] = generateHistogramOutline(theBins, theFrequency, theBinWidth);
	if (plotWithNegativePolarity)
		yy = -yy;
	end

	if (flipXY)
		theBarPlotHandle = RGCMosaicAnalyzer.visualize.shadedAreaBetweenTwoLines(ax, yy, xx, xx*0, ...
	        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);
	else
	    theBarPlotHandle = RGCMosaicAnalyzer.visualize.shadedAreaBetweenTwoLines(ax, xx, yy, yy*0, ...
	        faceColor, edgeColor, faceAlpha, lineWidth, lineStyle);
	end
end


function [xx, yy] = generateHistogramOutline(x,y, binWidth)

    xx(1) = x(1);
    yy(1) = 0;

    for i = 1:(numel(x))
    	
    	if (i < numel(x))
    		if (isempty(binWidth))
    			adaptiveBinWidth = (x(i+1)-x(i));
    			xOffset = 0.5*((x(i+1)-x(i)) - adaptiveBinWidth);
    		else
    			xOffset = 0.5*((x(i+1)-x(i)) - binWidth);
    		end
    	else
    		xOffset = 0.5*((x(i)-x(i-1)) - binWidth);
    	end

    	xx(numel(xx)+1) = x(i)+xOffset;
        yy(numel(yy)+1) = 0;

        xx(numel(xx)+1) = x(i)+xOffset;
        yy(numel(yy)+1) = y(i);

        if (i < numel(x))
	        xx(numel(xx)+1) = x(i+1)-xOffset;
	        yy(numel(yy)+1) = y(i);

	        xx(numel(xx)+1) = x(i+1)-xOffset;
	        yy(numel(yy)+1) = 0;
	    else
	    	xx(numel(xx)+1) = x(i)+(x(2)-x(1))-xOffset;
	        yy(numel(yy)+1) = y(i);

	        xx(numel(xx)+1) = x(i)+(x(2)-x(1))-xOffset;
	        yy(numel(yy)+1) = 0;
    	end
    end

    xx(numel(xx)+1) = x(end);
    yy(numel(yy)+1) = 0;

    xx = xx - (x(2)-x(1));
end
