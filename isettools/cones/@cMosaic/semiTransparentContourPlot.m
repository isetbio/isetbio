% Function to generate a semitransparent controur plot
function semiTransparentContourPlot(axesHandle, xSupport, ySupport, zData, zLevels, cmap, alpha, contourLineColor, varargin)
    
    p = inputParser;
    p.addParameter('lineWidth', 1.0, @isscalar);
    p.addParameter('edgeAlpha', 0.5, @isscalar);
    p.parse(varargin{:});
    lineWidth = p.Results.lineWidth;
    edgeAlpha = p.Results.edgeAlpha;

    
    % Compute contours at desired Z-level
    C = contourc(xSupport, ySupport, zData, zLevels);
    dataPoints = size(C,2);
    startPoint = 1;
    iContour = 1;
    while (startPoint < dataPoints)
        theLevel = C(1,startPoint);
        theNormalizedLevels(iContour) = (theLevel-min(zLevels))/(max(zLevels)-min(zLevels));
        theLevelVerticesNum = C(2,startPoint);
        points{iContour} = startPoint+(1:theLevelVerticesNum);
        contourLengths(iContour) = theLevelVerticesNum;
        iContour = iContour + 1;
        startPoint = startPoint + theLevelVerticesNum+1;
    end

    % Plot contours starting with the longest one, then the second longest,
    % etc.
    [~,idx] = sort(contourLengths, 'descend');
    hold(axesHandle, 'on');
    for iSortedContour = 1:numel(idx)
        theNormalizedLevel = theNormalizedLevels(idx(iSortedContour));
        x = C(1, points{idx(iSortedContour)});
        y = C(2, points{idx(iSortedContour)});
        v = [x(:) y(:)];
        f = 1:numel(x);
        lutIndex = 1+round(theNormalizedLevel*(size(cmap,1)-1));
        if (ischar(alpha))
            if (strcmp(alpha, 'matchZLevel'))
                theAlpha = max([0 0.5-0.7*theNormalizedLevel]);
            else
                error('Unknown alpha: ''%s''.', alpha);
            end
        else
            theAlpha = alpha;
        end

        patch('Faces', f, 'Vertices', v, ...
            'EdgeColor', [0 0 0], ...
            'FaceColor', cmap(lutIndex,:), ...
            'FaceAlpha', theAlpha, ...
            'EdgeAlpha', edgeAlpha, ...
            'EdgeColor', contourLineColor, ...
            'LineStyle', '-', 'LineWidth', lineWidth, ...
            'Parent', axesHandle);
    end

end