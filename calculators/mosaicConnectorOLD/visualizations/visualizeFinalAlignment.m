function visualizeFinalAlignment(conePositionsMicrons, RGCRFPositionsMicrons, X1, X2, Y1, Y2, desiredConesToRGCratios)
        hFig = figure(2); clf;
        set(hFig, 'Name', 'RGC mosaic alignment to cone mosaic');
        theAxesGrid = plotlab.axesGrid(hFig, ...
            'leftMargin', 0.06, ...
            'bottomMargin', 0.06, ...
            'rightMargin', 0.01, ...
            'topMargin', 0.05);

        scatter(theAxesGrid{1,1}, conePositionsMicrons(:,1), conePositionsMicrons(:,2), 'b'); hold on;
        scatter(theAxesGrid{1,1}, RGCRFPositionsMicrons(:,1), RGCRFPositionsMicrons(:,2), 300, 'g');

        conesInMosaicPatch = size(conePositionsMicrons,1);
        rgcsInMosaicPatch = size(RGCRFPositionsMicrons,1);
        
        plot(theAxesGrid{1,1},[X1; X2], ...
             [Y1; Y2], 'k-', 'LineWidth', 1.5);
        title(theAxesGrid{1,1},sprintf('cone-to-RGC ratio: %2.2f (desired), %2.2f (actual)', ...
            mean(desiredConesToRGCratios), conesInMosaicPatch/rgcsInMosaicPatch));
        drawnow;
end

