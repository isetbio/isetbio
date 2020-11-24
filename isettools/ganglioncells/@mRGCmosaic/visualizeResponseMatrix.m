function visualizeResponseMatrix(obj, independentVariable, responseMatrix)

    figSize = [1300 500];
    hFig = figure(2); clf;
    set(hFig, 'Position', [0 0 figSize(1) figSize(2)]);
    axesPosition = obj.axesMatrixPosition();
    for iRGC = 1:size(responseMatrix,1)
        ax = axes('Position', axesPosition(iRGC,:), 'Color', [1 1 1]);
        plot(ax,independentVariable, responseMatrix(iRGC,:), 'ko-');
        axis(ax, 'square');
        set(ax, 'XTick', [], 'YTick', [], 'XScale', 'log');
        drawnow
    end
    xlabel(ax,'spatial frequency (c/deg)');
    ylabel(ax,'response modulation');
    
end

