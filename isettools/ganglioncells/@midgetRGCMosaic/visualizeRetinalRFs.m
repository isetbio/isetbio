function visualizeRetinalRFs(obj, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('exportGraphicForEachRF', false, @islogical);
    p.addParameter('maxExportedGraphs', 15, @isnumeric);
    p.parse(varargin{:});

    exportGraphicForEachRF = p.Results.exportGraphicForEachRF;
    maxExportedGraphs = p.Results.maxExportedGraphs;


    mRGCmosaicCenterDegs = mean(obj.rgcRFpositionsDegs,1);
    mRGCsNum = size(obj.rgcRFcenterConeConnectivityMatrix,2);
    theRetinalRFmaps = cell(1, mRGCsNum);

    marginDegs = min([0.5 0.2*min(obj.sizeDegs)]);
    xLimsVisualized(1) = mRGCmosaicCenterDegs(1)-marginDegs;
    xLimsVisualized(2) = xLimsVisualized(1)+2*marginDegs;
    yLimsVisualized(1) = mRGCmosaicCenterDegs(2)-marginDegs;
    yLimsVisualized(2) = yLimsVisualized(1)+2*marginDegs;

   parfor iRGC = 1:mRGCsNum
        % Find the center cone indices
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        inputConeIndices = find(connectivityVector > 0.0001);
        allInputConeWeights = connectivityVector(inputConeIndices);
        allInputConePositionsDegs = obj.inputConeMosaic.coneRFpositionsDegs(inputConeIndices,:);
        allInputConeRcDegs = obj.inputConeMosaic.coneApertureDiametersDegs(inputConeIndices) * ...
                             obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;

        % Determine x and y range based on this mRGC center cone inputs
        minXY = min(allInputConePositionsDegs,[],1);
        maxXY = max(allInputConePositionsDegs,[],1);
        meanXY = mean(allInputConePositionsDegs,1);
        xRange = maxXY(1)-minXY(1);
        yRange = maxXY(2)-minXY(2);
        xyRange = max([xRange yRange]);
        x1 = meanXY(1)-xyRange/2-1.2*marginDegs;
        x2 = meanXY(1)+xyRange/2+1.2*marginDegs;
        y1 = meanXY(2)-xyRange/2-1.2*marginDegs;
        y2 = meanXY(2)+xyRange/2+1.2*marginDegs;

        spatialSupportDegsX = linspace(x1,x2, 200);
        spatialSupportDegsY = linspace(y1,y2, 200);
        [X,Y] = meshgrid(spatialSupportDegsX, spatialSupportDegsY);
        
        theRFcenterMap = X*0;
        for iInput = 1:numel(inputConeIndices)
             inputConePosDegs = allInputConePositionsDegs(iInput,:);
             inputConeRcDegs = allInputConeRcDegs(iInput);
             XX = X - inputConePosDegs(1);
             YY = Y - inputConePosDegs(2);
             R = sqrt(XX.^2 + YY.^2);
             theConeAperture = exp(-(R/inputConeRcDegs).^2);
             theRFcenterMap = theRFcenterMap + allInputConeWeights(iInput)*theConeAperture;
        end 

         theRetinalRFmaps{iRGC} = struct(...
             'centerRF', theRFcenterMap, ...
             'spatialSupportDegsX', spatialSupportDegsX, ...
             'spatialSupportDegsY', spatialSupportDegsY...
             );
    end % iRGC

    if (obj.sizeDegs < 0.4)
        tickIncrement = 0.05;
    elseif (obj.sizeDegs < 0.8)
        tickIncrement = 0.1;
    elseif (obj.sizeDegs < 1.6)
        tickIncrement = 0.2;
    elseif (obj.sizeDegs < 2)
        tickIncrement = 0.4;
    else
        tickIncrement = 0.5;
    end

    f = 1/tickIncrement;

    xMinTick = round(mRGCmosaicCenterDegs(1)*f)/f - 0.5*obj.sizeDegs(1);
    xMaxTick = round(mRGCmosaicCenterDegs(1)*f)/f + 0.5*obj.sizeDegs(1);
    xTicks = xMinTick:tickIncrement:xMaxTick;
    yMinTick = round(mRGCmosaicCenterDegs(2)*f)/f - 0.5*obj.sizeDegs(2);
    yMaxTick = round(mRGCmosaicCenterDegs(2)*f)/f + 0.5*obj.sizeDegs(2);
    yTicks = yMinTick:tickIncrement:yMaxTick;


    hFig = figure(1);
    set(hFig, 'Position', [10 10 880 880], 'Color', [1 1 1]);
    clf;

    % Sort RGCs
    ecc = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
    [~,sortedRGCindices] = sort(ecc, 'ascend');

    for iRGC = 1:numel(sortedRGCindices)

        if (exportGraphicForEachRF) && (iRGC > maxExportedGraphs)
            continue;
        end

        targetRGCindex  = sortedRGCindices(iRGC);
        %[D, idx] = MosaicConnector.pdist2(...
        %    obj.rgcRFpositionsDegs, obj.rgcRFpositionsDegs(targetRGCindex,:), 'smallest', 2);
        %theClosestNeighboringRGC = idx(2);

        s = theRetinalRFmaps{targetRGCindex};
        xLims = [s.spatialSupportDegsX(1) s.spatialSupportDegsX(end)];
        yLims = [s.spatialSupportDegsY(1) s.spatialSupportDegsY(end)];
        theRF = s.centerRF;
        theRFprofileX = sum(theRF,1);
        theRFprofileY = sum(theRF,2);
        maxProfile = max([max(theRFprofileX(:)) max(theRFprofileY(:))]);
        
        theRFprofileX  = theRFprofileX / maxProfile;
        theRFprofileY  = theRFprofileY / maxProfile;

        theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
            s.spatialSupportDegsX, s.spatialSupportDegsY, theRF, ...
            'flatTopGaussian', ~true, ...
            'forcedOrientationDegs', [], ...
            'rangeForEllipseRcYRcXratio', [1/1.4 1.4], ...
            'forcedCentroidXYpos', obj.rgcRFpositionsDegs(targetRGCindex,:), ...
            'globalSearch', true, ...
            'multiStartsNum', 8);
        fittedEllipsoidMap = theFittedGaussian.ellipsoidMap;
  
        fittedEllipsoidMapProfileX = sum(fittedEllipsoidMap,1);
        fittedEllipsoidMapProfileY = sum(fittedEllipsoidMap,2);
        maxFittedProfile = max([max(fittedEllipsoidMapProfileX(:)) max(fittedEllipsoidMapProfileY(:))]);

        fittedEllipsoidMapProfileX  = fittedEllipsoidMapProfileX/maxFittedProfile;
        fittedEllipsoidMapProfileY  = fittedEllipsoidMapProfileY/maxFittedProfile;
        

        zLevels = 0.05:0.1:0.9;

        ax = subplot(2,2,1);
        cla(ax);
        % Fitted ellipsoid at 1,2,3 sigma
        [~,c] = contourf(ax,s.spatialSupportDegsX, s.spatialSupportDegsY, fittedEllipsoidMap/max(fittedEllipsoidMap(:)), ...
            zLevels);
        c.LineWidth = 1.0;
        c.Color = 'none';
        hold(ax, 'on');

        % Positions of cones within this patch
        coneXpositionDegs = obj.inputConeMosaic.coneRFpositionsDegs(:,1);
        coneYpositionDegs = obj.inputConeMosaic.coneRFpositionsDegs(:,2);
        coneIndicesInNeigborhood = find(...
            (coneXpositionDegs>=xLims(1)) & ...
            (coneXpositionDegs<=xLims(2)) & ...
            (coneYpositionDegs>=yLims(1)) & ...
            (coneYpositionDegs<=yLims(2)) ...
            );

        
        for iCone = 1:numel(coneIndicesInNeigborhood)
            coneIndex = coneIndicesInNeigborhood(iCone);
            xo = obj.inputConeMosaic.coneRFpositionsDegs(coneIndex,1);
            yo = obj.inputConeMosaic.coneRFpositionsDegs(coneIndex,2);
            r = obj.inputConeMosaic.coneApertureDiametersDegs(coneIndex) * ...
                obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
            xx = xo + r * cosd(0:10:360);
            yy = yo + r * sind(0:10:360);
            switch (obj.inputConeMosaic.coneTypes(coneIndex))
                case cMosaic.LCONE_ID
                    faceColor = [1 0.2 0.5];
                case cMosaic.MCONE_ID
                    faceColor = [0.2 1 0.5];
                case cMosaic.SCONE_ID
                    faceColor = [0.5 0.1 0.9];
            end

            patch(xx, yy, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.4, ...
                'FaceColor', faceColor, 'LineWidth', 0.5, ...
                'LineStyle', '-', ...
                'Parent', ax);
        end

        % RF map
        [~,c] = contour(ax,s.spatialSupportDegsX, s.spatialSupportDegsY, theRF, zLevels);
        c.LineWidth = 1.0;
        c.LineColor = 'black';

        % Cross-hairs
        plot(ax, xLimsVisualized, obj.rgcRFpositionsDegs(targetRGCindex,2)*[1 1], 'b-', 'LineWidth', 1.0);
        plot(ax, obj.rgcRFpositionsDegs(targetRGCindex,1)*[1 1], yLims, 'b-', 'LineWidth', 1.0);
        hold(ax, 'off');

        % Finalize graph
        set(ax, 'CLim', [0 1]);
        axis(ax, 'equal');
        set(ax, 'XLim', xLimsVisualized, 'YLim', yLimsVisualized);
        set(ax, 'XTick', xTicks, 'YTick', yTicks);
        set(ax, 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5],  'Color', [0 0 0], 'FontSize', 16);
        colormap(ax,brewermap(1024, '*greys'));

        ax = subplot(2,2,4);
        [~,c] = contourf(ax,s.spatialSupportDegsX, s.spatialSupportDegsY, fittedEllipsoidMap/max(fittedEllipsoidMap(:)), zLevels);
        c.LineWidth = 1.0;
        c.Color = 'none';
        axis(ax, 'equal');
        set(ax, 'XLim', mRGCmosaicCenterDegs(1) + 0.5*obj.sizeDegs(1) * [-1 1], ...
                'YLim', mRGCmosaicCenterDegs(2) + 0.5*obj.sizeDegs(2) * [-1 1], ...
                'XTick', xTicks, 'YTick', yTicks, 'CLim', [0 1]);
        set(ax, 'GridColor', [1 1 1],  'Color', [0 0 0], 'FontSize', 16, 'LineWidth', 1.0);
        grid(ax, 'on');
        colormap(ax,brewermap(1024, '*greys'));




        x = s.spatialSupportDegsX;
        y = s.spatialSupportDegsY;
        [X,Y] = meshgrid(x,y);
        
        inputConeIndices = find(obj.rgcRFcenterConeConnectivityMatrix(:,targetRGCindex)>0.0001);
        coneApertures = zeros(numel(inputConeIndices), size(X,1), size(X,2));
        for iCone = 1:numel(inputConeIndices)
            coneIndex = inputConeIndices(iCone);
            coneWeight = full(obj.rgcRFcenterConeConnectivityMatrix(coneIndex,targetRGCindex));

            xo = obj.inputConeMosaic.coneRFpositionsDegs(coneIndex,1);
            yo = obj.inputConeMosaic.coneRFpositionsDegs(coneIndex,2);
            rc = obj.inputConeMosaic.coneApertureDiametersDegs(coneIndex) * ...
                obj.inputConeMosaic.coneApertureToConeCharacteristicRadiusConversionFactor;
            
            R = sqrt((X-xo).^2+(Y-yo).^2);
            coneApertures(iCone,:,:) = coneWeight * exp(-(R/rc).^2);
        end

        % Y-profile
        ax = subplot(2,2,2);
        cla(ax);
        
        patch(fittedEllipsoidMapProfileY, s.spatialSupportDegsY, [0 0 0], ...
                'EdgeColor', [0 0 0], 'FaceAlpha', 0.4, ...
                'FaceColor', [0.8 0.8 0.8], 'LineWidth', 1, ...
                'LineStyle', '--', ...
                'Parent', ax);
        hold(ax, 'on');

        % Gaussian cone apertures
        for iCone = 1:numel(inputConeIndices)
            coneIndex = inputConeIndices(iCone);
            coneAperture = squeeze(coneApertures(iCone,:,:));
            coneApertureProfile = sum(coneAperture,2)/maxProfile;
            switch (obj.inputConeMosaic.coneTypes(coneIndex))
                case cMosaic.LCONE_ID
                    faceColor = [1 0.2 0.5];
                case cMosaic.MCONE_ID
                    faceColor = [0.2 1 0.5];
                case cMosaic.SCONE_ID
                    faceColor = [0.5 0.1 0.9];
            end
            shadedAreaPlot(ax,coneApertureProfile,y, 0, faceColor, faceColor*0.5, 0.5, 1.0, '-');
            
        end

        plot(ax, theRFprofileY, s.spatialSupportDegsY, 'k-', 'LineWidth', 1.5);
        
        % Cross-hairs
        plot(ax, [0 1], obj.rgcRFpositionsDegs(targetRGCindex,2)*[1 1], 'b-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax, 'square'); grid(ax, 'off');box(ax, 'off');
        set(ax, 'XColor', 'none', 'YColor', 'none', 'YLim', yLimsVisualized, 'YTick', yTicks, 'XLim', [0 1], 'XTick', 0:0.1:1);

        % X-profile
        ax = subplot(2,2,3);
        cla(ax);

        patch(s.spatialSupportDegsX, fittedEllipsoidMapProfileX, [0 0 0], ...
                'EdgeColor', [0 0 0], 'FaceAlpha', 0.4, ...
                'FaceColor', [0.8 0.8 0.8], 'LineWidth', 1, ...
                'LineStyle', '--', ...
                'Parent', ax);

        hold(ax, 'on')
        % Gaussian cone apertures
        for iCone = 1:numel(inputConeIndices)
            coneIndex = inputConeIndices(iCone);
            coneAperture = squeeze(coneApertures(iCone,:,:));
            coneApertureProfile = sum(coneAperture,1)/maxProfile;
            switch (obj.inputConeMosaic.coneTypes(coneIndex))
                case cMosaic.LCONE_ID
                    faceColor = [1 0.2 0.5];
                case cMosaic.MCONE_ID
                    faceColor = [0.2 1 0.5];
                case cMosaic.SCONE_ID
                    faceColor = [0.5 0.1 0.9];
            end
            shadedAreaPlot(ax,x,coneApertureProfile, 0, faceColor, faceColor*0.5, 0.5, 1.0, '-');
        end
        

        plot(ax, s.spatialSupportDegsX, theRFprofileX, 'k-', 'LineWidth', 1.5);
        % Cross-hairs
        plot(ax, obj.rgcRFpositionsDegs(targetRGCindex,1)*[1 1], [0 1], 'b-', 'LineWidth', 1.0);
        hold(ax, 'off');
        axis(ax, 'square'); grid(ax, 'off'); box(ax, 'off');
        set(ax, 'XColor', 'none', 'YColor', 'none', 'XLim', xLimsVisualized, 'XTick', xTicks, 'YLim', [0 1], 'YTick', 0:0.1:1);
        
        drawnow;
        if (exportGraphicForEachRF) && (iRGC <= maxExportedGraphs)
            fName = sprintf('RGC_%d_overlap_%2.2f.pdf', iRGC, obj.rfOverlapRatio);
            NicePlot.exportFigToPDF(fName, hFig, 300);
        end
    end

end



function shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)

    x = reshape(x, [1 numel(x)]);
    y = reshape(y, [1 numel(y)]);
    x = [x(1)      x  fliplr(x)     x(1) ];
    y = [baseline  y  y*0+baseline  baseline];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end