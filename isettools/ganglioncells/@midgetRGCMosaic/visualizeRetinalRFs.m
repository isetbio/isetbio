function visualizeRetinalRFs(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('exportGraphicForEachRF', false, @islogical);
    p.addParameter('maxExportedGraphs', 15, @isnumeric);
    p.addParameter('spatialSupportSamplesNum', 256, @isnumeric);
    p.parse(varargin{:});
    exportGraphicForEachRF = p.Results.exportGraphicForEachRF;
    maxExportedGraphs = p.Results.maxExportedGraphs;
    spatialSupportSamplesNum = p.Results.spatialSupportSamplesNum;

    % Compute additional margin for the spatial support
    mRGCmosaicCenterDegs = mean(obj.rgcRFpositionsDegs,1);
    marginDegs = min([0.5 0.4*min(obj.sizeDegs)]);

    % Compute visualization limints
    xLimsVisualized(1) = mRGCmosaicCenterDegs(1)-marginDegs;
    xLimsVisualized(2) = xLimsVisualized(1)+2*marginDegs;
    yLimsVisualized(1) = mRGCmosaicCenterDegs(2)-marginDegs;
    yLimsVisualized(2) = yLimsVisualized(1)+2*marginDegs;

    minXY = min(obj.rgcRFpositionsDegs,[],1);
    maxXY = max(obj.rgcRFpositionsDegs,[],1);
    xLimsVisualizedFull(1) = minXY(1)-max(obj.rgcRFspacingsDegs);
    xLimsVisualizedFull(2) = maxXY(1)+max(obj.rgcRFspacingsDegs);
    yLimsVisualizedFull(1) = minXY(2)-max(obj.rgcRFspacingsDegs);
    yLimsVisualizedFull(2) = maxXY(2)+max(obj.rgcRFspacingsDegs);

    xLimsVisualizedFull
    yLimsVisualizedFull
    % Compute the retinal RFcenter maps
    theRetinalRFcenterMaps = obj.computeRetinalRFcenterMaps(marginDegs, spatialSupportSamplesNum);

    % Sort RGCs according to their eccentricity
    ecc = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
    [~,sortedRGCindices] = sort(ecc, 'ascend');

     % Format figure
    hFig = figure(1);
    set(hFig, 'Position', [10 10 880 880], 'Color', [0.8 0.8 0.8]);
    
    hFig2 = figure(2); clf;
    set(hFig2, 'Position', [10 10 880 880], 'Color', [1 1 1]);
    ax2 = subplot('Position', [0.09 0.09 0.89 0.89]);
    obj.inputConeMosaic.visualize(...
        'figureHandle', hFig2, ...
        'axesHandle', ax2, ...
        'domain', 'degrees', ...
        'visualizedConeAperture', 'geometricArea', ...
        'visualizedConeApertureThetaSamples', 30, ...
        'backgroundColor', [1 1 1]);
    hold(ax2, 'on');
    drawnow;

    % Fitting the discrete RF center cone map
    fitTheDiscreteRFcenterMap = true;

    % Plot the RF maps
    for iRGC = 1:numel(sortedRGCindices)

        % Retrieve the RGCindex
        targetRGCindex  = sortedRGCindices(iRGC);

        % Retrieve the computed retinal center RF map
        s = theRetinalRFcenterMaps{targetRGCindex};
        theRF = s.centerRF;
        xLims = [s.spatialSupportDegsX(1) s.spatialSupportDegsX(end)];
        yLims = [s.spatialSupportDegsY(1) s.spatialSupportDegsY(end)];

        % Compute X and Y spatial profiles
        theRFprofileX = sum(theRF,1);
        theRFprofileY = sum(theRF,2);
        maxProfile = max([max(theRFprofileX(:)) max(theRFprofileY(:))]);
        theRFprofileX  = theRFprofileX / maxProfile;
        theRFprofileY  = theRFprofileY / maxProfile;

        if (fitTheDiscreteRFcenterMap)
            % Fit the discrete center RF map with an ellipsoidal Gaussian
            theFittedGaussian = RetinaToVisualFieldTransformer.fitScatterGaussianEllipsoid(...
                s.spatialSupportDegsX, s.spatialSupportDegsY, theRF,...
                s.inputConeWeights, obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices,:), ...
                'flatTopGaussian', ~true, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/1.4 1.4], ...
                'forcedCentroidXYpos', obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                'globalSearch', true, ...
                'multiStartsNum', 8);
        else
            % Fit the continuous center RF map with an ellipsoidal Gaussian
            theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                s.spatialSupportDegsX, s.spatialSupportDegsY, theRF, ...
                'flatTopGaussian', ~true, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/1.4 1.4], ...
                'forcedCentroidXYpos', obj.rgcRFpositionsDegs(targetRGCindex,:), ...
                'globalSearch', true, ...
                'multiStartsNum', 8);
        end

        % Extract the fitted ellipsoidal Gaussian
        fittedEllipsoidMap = theFittedGaussian.ellipsoidMap;

        % Compute X and Y spatial profiles of the fitted ellipsoidal Gaussian
        fittedEllipsoidMapProfileX = sum(fittedEllipsoidMap,1);
        fittedEllipsoidMapProfileY = sum(fittedEllipsoidMap,2);
        maxFittedProfile = max([max(fittedEllipsoidMapProfileX(:)) max(fittedEllipsoidMapProfileY(:))]);
        fittedEllipsoidMapProfileX  = fittedEllipsoidMapProfileX/maxFittedProfile;
        fittedEllipsoidMapProfileY  = fittedEllipsoidMapProfileY/maxFittedProfile;
        
        
        if (iRGC <= maxExportedGraphs)
            % Clear the figure
            figure(hFig);
            clf;
    
            % The RF map together with the fitted ellipsoid Gaussian and cones
            ax = subplot(2,2,1);
    
            % Contour of the fitted ellipsoid
            zLevels = 0.05:0.1:0.9;
            [~,c] = contourf(ax,s.spatialSupportDegsX, s.spatialSupportDegsY, fittedEllipsoidMap/max(fittedEllipsoidMap(:)), zLevels);
            c.LineWidth = 1.0; c.Color = 'none';
            hold(ax, 'on');
    
            % Find the indices of cones in the neighborhood
            plotConesInTheNeighborhood(ax, obj.inputConeMosaic, xLims, yLims)
    
            % Compute ticks
            if (iRGC == 1)
                [xTicks, yTicks] = computeTicks(obj.sizeDegs, mRGCmosaicCenterDegs);
            end
    
            % Contour of the RF map
            [~,c] = contour(ax,s.spatialSupportDegsX, s.spatialSupportDegsY, theRF, zLevels);
            c.LineWidth = 1.0;
            c.LineColor = 'black';
    
            % Cross-hairs in blue
            plot(ax, xLimsVisualized, obj.rgcRFpositionsDegs(targetRGCindex,2)*[1 1], 'c-', 'LineWidth', 1.0);
            plot(ax, obj.rgcRFpositionsDegs(targetRGCindex,1)*[1 1], yLims, 'c-', 'LineWidth', 1.0);
            hold(ax, 'off');
    
            % Finalize graph
            set(ax, 'CLim', [0 1]);
            axis(ax, 'equal');
            set(ax, 'XLim', xLimsVisualized, 'YLim', yLimsVisualized);
            set(ax, 'XTick', xTicks, 'YTick', yTicks);
            set(ax, 'XColor', [0.25 0.25 0.25], 'YColor', [0.25 0.25 0.25],  'Color', [0 0 0], 'FontSize', 16, 'LineWidth', 1.0);
            xlabel('eccentricity, x (degs)')
            ylabel('eccentricity, y (degs)')
            colormap(ax,brewermap(1024, '*greys'));
            
    
            % The fitted ellipsoidal Gaussian alone
            ax = subplot(2,2,4);
            [~,c] = contourf(ax,s.spatialSupportDegsX, s.spatialSupportDegsY, fittedEllipsoidMap/max(fittedEllipsoidMap(:)), zLevels);
            c.LineWidth = 1.0; c.Color = 'none';
            axis(ax, 'equal');
            set(ax, 'XLim', mRGCmosaicCenterDegs(1) + 0.5*obj.sizeDegs(1) * [-1 1], ...
                    'YLim', mRGCmosaicCenterDegs(2) + 0.5*obj.sizeDegs(2) * [-1 1], ...
                    'XTick', xTicks, 'YTick', yTicks, 'CLim', [0 1]);
            set(ax, 'XColor', [0.25 0.25 0.25], 'YColor', [0.25 0.25 0.25],  'Color', [0 0 0], 'FontSize', 16, 'LineWidth', 1.0);
            grid(ax, 'on');
            xlabel(ax,'eccentricity, x (degs)')
            ylabel(ax,'eccentricity, y (degs)')
            colormap(ax,brewermap(1024, '*greys'));
    
    
    
            % The Y-profile of the fitted Gaussian ellipsoid
            ax = subplot(2,2,2);
    
            % The fitted ellipsoid Y-profile in gray
            patch(fittedEllipsoidMapProfileY, s.spatialSupportDegsY, [0 0 0], ...
                    'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 1.0, ...
                    'FaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5, ...
                    'LineStyle', '-', ...
                    'Parent', ax);
            hold(ax, 'on');
    
            % The input cone aperture Y-profiles
            for iCone = 1:numel(s.inputConeIndices)
                coneIndex = s.inputConeIndices(iCone);
                coneAperture = squeeze(s.coneApertures(iCone,:,:));
                coneApertureProfile = sum(coneAperture,2)/maxProfile;
                faceColor = coneColor(obj.inputConeMosaic.coneTypes(coneIndex));
                shadedAreaPlot(ax,coneApertureProfile,s.spatialSupportDegsY, 0, faceColor, faceColor*0.5, 0.5, 1.0, '-');
            end
    
            % The actual Y-profile of the retinal RF
            plot(ax, theRFprofileY, s.spatialSupportDegsY, 'y-', 'LineWidth', 1.5);
            
            % Cross-hairs
            plot(ax, [0 1], obj.rgcRFpositionsDegs(targetRGCindex,2)*[1 1], 'c-', 'LineWidth', 1.0);
            hold(ax, 'off');
            axis(ax, 'square'); grid(ax, 'off');box(ax, 'off');
            set(ax, 'XColor', 'none', 'YColor', 'none', 'YLim', yLimsVisualized, 'YTick', yTicks, 'XLim', [0 1], 'XTick', 0:0.1:1);
            set(ax, 'Color', [0 0 0], 'FontSize', 16, 'LineWidth', 1.0);
            ylabel(ax,'eccentricity, y (degs)')
    
    
            % The X-profile of the fitted Gaussian ellipsoid
            ax = subplot(2,2,3);
     
            % The fitted ellipsoid X- profile in gray
            patch(s.spatialSupportDegsX, fittedEllipsoidMapProfileX, [0 0 0], ...
                    'EdgeColor', [0.3 0.3 0.3], 'FaceAlpha', 1.0, ...
                    'FaceColor', [0.8 0.8 0.8], 'LineWidth', 1.5, ...
                    'LineStyle', '--', ...
                    'Parent', ax);
    
            hold(ax, 'on')
    
            % The input cone aperture X-profiles
            for iCone = 1:numel(s.inputConeIndices)
                coneIndex = s.inputConeIndices(iCone);
                coneAperture = squeeze(s.coneApertures(iCone,:,:));
                coneApertureProfile = sum(coneAperture,1)/maxProfile;
                faceColor = coneColor(obj.inputConeMosaic.coneTypes(coneIndex));
                shadedAreaPlot(ax,s.spatialSupportDegsX, coneApertureProfile, 0, faceColor, faceColor*0.5, 0.5, 1.0, '-');
            end
            
            % The actual X-profile of the retinal RF
            plot(ax, s.spatialSupportDegsX, theRFprofileX, 'y-', 'LineWidth', 1.5);
    
            % Cross-hairs
            plot(ax, obj.rgcRFpositionsDegs(targetRGCindex,1)*[1 1], [0 1], 'c-', 'LineWidth', 1.0);
            hold(ax, 'off');
            axis(ax, 'square'); grid(ax, 'off'); box(ax, 'off');
            set(ax, 'XColor', 'none', 'YColor', 'none', 'XLim', xLimsVisualized, 'XTick', xTicks, 'YLim', [0 1], 'YTick', 0:0.1:1);
            set(ax, 'Color', [0 0 0], 'FontSize', 16, 'LineWidth', 1.0);
            xlabel(ax,'eccentricity, x (degs)')
    
            % Export to PDF
            drawnow;
            if (exportGraphicForEachRF)
                fName = sprintf('RGC_%d_Overlap_%2.2f_Ecc_%2.1f_%2.1f.pdf', iRGC, obj.rfOverlapRatio, mRGCmosaicCenterDegs(1), mRGCmosaicCenterDegs(2));
                NicePlot.exportFigToPDF(fName, hFig, 300);
            end
        end

        
        figure(hFig2);
        oneSigmaLevel = [1 exp(-0.5)];
        cMap = brewermap(128, 'greys');
        alpha = 0.5;
        contourLineColor = [0 0 0];
        cMosaic.semiTransparentContourPlot(ax2, s.spatialSupportDegsX, s.spatialSupportDegsY, fittedEllipsoidMap/max(fittedEllipsoidMap(:)), ...
            oneSigmaLevel, cMap, alpha, contourLineColor, ...
            'lineWidth', 2.0, ...
            'edgeAlpha', 1.0);
    
        hold(ax2, 'on');


        % Lines connecting the centroid with all the inputs with line
        % darkness indicating weight: the stronger the weight, the darker the line
        if (numel(s.inputConeIndices)>1)
            inputsCentroid = mean(obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices,:),1);
            for iInput = 1:numel(s.inputConeIndices)
                if (s.inputConeWeights(iInput)>0.001)
                    plot(ax2,[inputsCentroid(1) obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices(iInput),1)], ...
                        [inputsCentroid(2) obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices(iInput),2)], ...
                        'k-', 'LineWidth', 1.0, 'Color', [1 1 1]*(1-s.inputConeWeights(iInput)/max(s.inputConeWeights(:))));
                end
            end
        end

        
        axis(ax2, 'equal');
        set(ax2, 'XLim', xLimsVisualizedFull, ...
                 'YLim', yLimsVisualizedFull, ...
                 'XTick', xTicks, 'YTick', yTicks, ...
                 'XTickLabel', sprintf('%2.1f\n', xTicks), 'YTickLabel', sprintf('%2.1f\n', yTicks), ...
                 'CLim', [0 1]);
        set(ax2, 'XColor', [0.25 0.25 0.25], 'YColor', [0.25 0.25 0.25], 'FontSize', 20, 'LineWidth', 1.0);
        grid(ax2, 'on');
        xlabel(ax2,'eccentricity, x (degs)')
        ylabel(ax2,'eccentricity, y (degs)')
        drawnow;
        

    end

    if (exportGraphicForEachRF)
        fName = sprintf('allRGCsOneSigmaOutline_Overlap_%2.2f_Ecc_%2.1f_%2.1f.pdf', obj.rfOverlapRatio, mRGCmosaicCenterDegs(1), mRGCmosaicCenterDegs(2));
        NicePlot.exportFigToPDF(fName, hFig2, 300);
    end


end


% Supporting functions
function theColor = coneColor(coneType)
    switch (coneType)
        case cMosaic.LCONE_ID
            theColor = [1 0.2 0.5];
        case cMosaic.MCONE_ID
            theColor = [0.2 1 0.5];
        case cMosaic.SCONE_ID
            theColor = [0.5 0.1 0.9];
    end
end


function plotConesInTheNeighborhood(ax, inputConeMosaic, xLims, yLims)

    coneXpositionDegs = inputConeMosaic.coneRFpositionsDegs(:,1);
    coneYpositionDegs = inputConeMosaic.coneRFpositionsDegs(:,2);
    coneIndicesInNeigborhood = find(...
        (coneXpositionDegs>=xLims(1)) & ...
        (coneXpositionDegs<=xLims(2)) & ...
        (coneYpositionDegs>=yLims(1)) & ...
        (coneYpositionDegs<=yLims(2)) ...
        );

    % Plot all cones in the neighborhood
    for iCone = 1:numel(coneIndicesInNeigborhood)
        coneIndex = coneIndicesInNeigborhood(iCone);
        xo = inputConeMosaic.coneRFpositionsDegs(coneIndex,1);
        yo = inputConeMosaic.coneRFpositionsDegs(coneIndex,2);
        r = inputConeMosaic.coneApertureDiametersDegs(coneIndex) * ...
            inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
        xx = xo + r * cosd(0:10:360);
        yy = yo + r * sind(0:10:360);
        switch (inputConeMosaic.coneTypes(coneIndex))
            case cMosaic.LCONE_ID
                faceColor = [1 0.2 0.5];
            case cMosaic.MCONE_ID
                faceColor = [0.2 1 0.5];
            case cMosaic.SCONE_ID
                faceColor = [0.5 0.1 0.9];
        end

        patch(xx, yy, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
            'FaceColor', faceColor, 'LineWidth', 0.5, ...
            'LineStyle', '-', ...
            'Parent', ax);
    end
end


function  [xTicks, yTicks] = computeTicks(mosaicSize, mosaicCenter)
    if (mosaicSize< 0.4)
        tickIncrement = 0.05;
    elseif (mosaicSize < 0.8)
        tickIncrement = 0.1;
    elseif (mosaicSize < 1.6)
        tickIncrement = 0.2;
    elseif (mosaicSize < 2)
        tickIncrement = 0.4;
    else
        tickIncrement = 0.5;
    end

    f = 1/tickIncrement;

    xMinTick = round(mosaicCenter(1)*f)/f - 0.5*mosaicSize(1);
    xMaxTick = round(mosaicCenter(1)*f)/f + 0.5*mosaicSize(1);
    xTicks = xMinTick:tickIncrement:xMaxTick;
    yMinTick = round(mosaicCenter(2)*f)/f - 0.5*mosaicSize(2);
    yMaxTick = round(mosaicCenter(2)*f)/f + 0.5*mosaicSize(2);
    yTicks = yMinTick:tickIncrement:yMaxTick;
end



function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = reshape(x, [1 numel(x)]);
    y = reshape(y, [1 numel(y)]);
    px = [x(1)      x  fliplr(x)     x(1) ];
    py = [baseline  y  y*0+baseline  baseline];
    ii = find(abs(py)> 100*eps);
    px = px(ii);
    py = py(ii);
    pz = -1*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end