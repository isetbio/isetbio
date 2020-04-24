function plotMosaic(theAxesHandle, rfPositions, visualizedFOVMicrons, mosaicTitle)
    if (~isempty(rfPositions))
        plot(theAxesHandle, rfPositions(:,1), rfPositions(:,2), 'k.'); 
        axis(theAxesHandle, 'equal');
        xRange = 0.5*visualizedFOVMicrons*[-1 1] + 0.5*(max(squeeze(rfPositions(:,1)))+min(squeeze(rfPositions(:,1))));
        yRange = 0.5*visualizedFOVMicrons*[-1 1] + 0.5*(max(squeeze(rfPositions(:,2)))+min(squeeze(rfPositions(:,2))));
        axis(theAxesHandle, 'square');
        set(theAxesHandle , 'XLim', xRange, 'YLim', yRange);
        xlabel(theAxesHandle, 'space (microns)');
        if (~isempty(mosaicTitle))
            title(theAxesHandle,mosaicTitle, 'Color', 'r');
        end
    end
end