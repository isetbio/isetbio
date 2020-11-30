function visualizeResponseMatrix(obj, independentVariable, responseMatrix, varargin)
    p = inputParser;
    p.addParameter('fittedResponses', @(x)(isempty(x) || isstruct(x)));
    p.parse(varargin{:}); 
    fittedResponses = p.Results.fittedResponses;
    
    figSize = [1300 1300];
    hFig = figure(2); clf;
    set(hFig, 'Position', [0 0 figSize(1) figSize(2)]);
    axesPosition = obj.axesMatrixPosition();
    for iRGC = 1:size(responseMatrix,1)
        ax = axes('Position', axesPosition(iRGC,:), 'Color', [1 1 1]);
        hold(ax, 'on');
        if (~isempty(fittedResponses))
            plot(ax,independentVariable, responseMatrix(iRGC,:), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5]);
            plot(fittedResponses.x, fittedResponses.y(iRGC,:), 'r-', 'Color', [0.6 0 0], 'LineWidth', 1.5);
        else
            plot(ax,independentVariable, responseMatrix(iRGC,:), 'ro-', 'MarkerSize', 8, 'MarkerFaceColor', [1 0.5 0.5]);
        end
        
        axis(ax, 'square');
        grid(ax, 'off');
        set(ax, 'XTick', [0.3 1 3 10 30], 'YTick', [], 'XScale', 'log');
        drawnow
    end
    xlabel(ax,'spatial frequency (c/deg)');
    ylabel(ax,'response modulation');
    
end

