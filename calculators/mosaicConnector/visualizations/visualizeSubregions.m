function visualizeSubregions(figNo,midgetRGCconnectionMatrixCenter, midgetRGCconnectionMatrixSurround, rgcEccDegs, ...
    rgcCenterPositionMicrons, retinalSurroundRadiiMicrons, conePositionsMicrons, coneSpacingsMicrons,  coneTypes, plotlabOBJ, pdfFileName,exportsDir)

    maxGainDivisor = 500;

    hFig = figure(figNo); clf;
    theAxesGrid = plotlabOBJ.axesGrid(hFig, ...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'leftMargin', 0.08, ...
        'widthMargin', 0.08, ...
        'heightMargin', 0.02, ...
        'bottomMargin', 0.06, ...
        'rightMargin', 0.01, ...
        'topMargin', 0.01);
    
   
    videoOBJ = VideoWriter(fullfile(exportsDir, pdfFileName), 'MPEG-4'); % H264 format
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();
        
    
    previousEccDegs = [];
    previousCenterRadii = [];
    previousSurroundRadii = [];
    previousIntegratedSensitivityRatios = [];
    
    rgcsNum = size(rgcCenterPositionMicrons,1);
    
    randomizedK = randperm(rgcsNum);
    for k = 1:numel(randomizedK)   
        
        % Get randomized index
        iRGC = randomizedK(k);
        
        % Get position of RGC in microns
        rgcPositionMicrons = rgcCenterPositionMicrons(iRGC,:);
         
        % center weights
        connectivityVector = full(squeeze(midgetRGCconnectionMatrixCenter(:, iRGC)));
        coneIndicesConnectedToCenter = find(connectivityVector>0);
        centerWeights = squeeze(full(midgetRGCconnectionMatrixCenter(coneIndicesConnectedToCenter, iRGC)));
        
        % surround weights
        connectivityVector = full(squeeze(midgetRGCconnectionMatrixSurround(:, iRGC)));
        coneIndicesConnectedToSurround = find(connectivityVector>0);
        surroundWeights = squeeze(full(midgetRGCconnectionMatrixSurround(coneIndicesConnectedToSurround, iRGC)));
        
        
        % Compute visualized spatial support
        rgcSurroundRadiusMicrons = retinalSurroundRadiiMicrons(iRGC);
        % Varied support, scaling with surround
        spatialSupportRangeMicrons = rgcSurroundRadiusMicrons*4; 
        % or fixed support
        %spatialSupportRangeMicrons = round(0.4*300);
        
        % Plot center weights
        % Only visualize weights up to a fraction of the mean(center). After that saturate to max
        maxWeightVisualized = max(centerWeights);
        
        visualizeWeigtedConeInputsToRGCSubregion(theAxesGrid{1,1}, coneTypes(coneIndicesConnectedToCenter), ...
            conePositionsMicrons(coneIndicesConnectedToCenter,:), ...
            coneSpacingsMicrons(coneIndicesConnectedToCenter), ...
            centerWeights, maxWeightVisualized, ...
            rgcPositionMicrons,  spatialSupportRangeMicrons, 1);
        
        % Plot surround weights
        % Only visualize weights up to a fraction of the mean(center). After that saturate to max
        maxWeightVisualized = max(surroundWeights);
        visualizeWeigtedConeInputsToRGCSubregion(theAxesGrid{1,2}, coneTypes(coneIndicesConnectedToSurround), ...
            conePositionsMicrons(coneIndicesConnectedToSurround,:), ...
            coneSpacingsMicrons(coneIndicesConnectedToSurround), ...
            surroundWeights,  maxWeightVisualized, ...
            rgcPositionMicrons, spatialSupportRangeMicrons, -1);
        
        if (1==2)
            % Plot the center&surround radii
            theAxes = theAxesGrid{2,1};
            cla(theAxes);
            hold(theAxes, 'on');
            if (~isempty(previousEccDegs))
                 scatter(theAxes, previousEccDegs, previousCenterRadii, 'o');
                 scatter(theAxes, previousEccDegs, previousSurroundRadii, 'd');
            end
        
           % scatter(theAxes, rgcEccDegs(iRGC), rgcCenterRadiusDegs, 225,'o', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [1 0 0]);
           % scatter(theAxes, rgcEccDegs(iRGC),  rgcSurroundRadiusDegs, 225,'d', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [0 0 1]);

            set(theAxes, 'XScale', 'log', 'XLim', [0.006 30], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30]);
            set(theAxes, 'YScale', 'log', 'YLim', [0.002 1], ...
                'YTick', [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10], 'YTickLabel', {'0.001', '0.003', '0.01', '0.03',' 0.1', '0.3', '1', '3', '10'});
            xlabel(theAxes,'eccentricity (degs)');
            ylabel(theAxes,'subregion radius (degs)');

            % Plot integrated sensitivity ratio
            theAxes = theAxesGrid{2,2};
            cla(theAxes);
            hold(theAxes, 'on');
            if (~isempty(previousEccDegs))
                scatter(theAxes, previousEccDegs, previousIntegratedSensitivityRatios, 'o', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0 0 0]);
            end
            scatter(theAxes, rgcEccDegs(iRGC), sum(surroundWeights)/sum(centerWeights), 225, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerFaceAlpha', 1.0, 'MarkerEdgeColor', [0 0 0]);

            set(theAxes, 'XScale', 'log', 'XLim', [0.006 30], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30], 'YLim', [0 1]);
            xlabel(theAxes, 'eccentricity (degs)');
            ylabel(theAxes, 'integrated sensitivity (surround/center)');

            previousEccDegs = cat(2, previousEccDegs, rgcEccDegs(iRGC));
            previousCenterRadii = cat(2, previousCenterRadii, rgcCenterRadiusDegs);
            previousSurroundRadii = cat(2, previousSurroundRadii, rgcSurroundRadiusDegs);
            previousIntegratedSensitivityRatios = cat(2, previousIntegratedSensitivityRatios, sum(surroundWeights)/sum(centerWeights));
        end
        
        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
         
        %plotlabOBJ.exportFig(hFig, 'png', pdfFileName,exportsDir);
    end
    videoOBJ.close();
end
