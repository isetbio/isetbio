function offsetAxes(ax, ff, XLims, YLims, varargin)
    
    p = inputParser;
    p.addParameter('keepXaxisSymmetric', false, @islogical);
    p.addParameter('keepYaxisSymmetric', false, @islogical);
    p.parse(varargin{:});
    keepXaxisSymmetric = p.Results.keepXaxisSymmetric;
    keepYaxisSymmetric = p.Results.keepYaxisSymmetric;

    figureWidthHeightRatio = ax.Parent.Position(3)/ax.Parent.Position(4);

    if (isempty(XLims))
        XLims = get(ax, 'XLim');
    end

    if (isempty(YLims))
        YLims = get(ax, 'YLim');
    end

    dX = XLims(2)-XLims(1);
    dY = YLims(2)-YLims(1);
    xOffset = dX*ff.axisOffsetFactor;
    yOffset = dY*ff.axisOffsetFactor;

    if (keepXaxisSymmetric)
        XLims = [XLims(1)+xOffset XLims(2)-xOffset];
    else
        XLims = [XLims(1)+xOffset XLims(2)];
    end

    if (keepYaxisSymmetric)
        YLims = [YLims(1)+yOffset YLims(2)-yOffset];
    else
        YLims = [YLims(1)+yOffset YLims(2)];
    end

    set(ax, 'XLim', XLims, 'YLim', YLims);
end