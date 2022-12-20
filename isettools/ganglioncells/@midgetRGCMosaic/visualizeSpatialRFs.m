function [hFig, allAxes] = visualizeSpatialRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('maxVisualizedRFs', 18, @(x)(isempty(x) || isscalar(x)));
    p.addParameter('onlyForRGCwithIndex', [], @(x)(isempty(x) || isscalar(x)));
    p.addParameter('generateVideo', false, @islogical);
    p.addParameter('videoFileName', [], @(x)(isempty(x) || ischar(x)));
    p.addParameter('withEccentricityCrossHairs', false, @islogical);
    p.addParameter('withPSFData', [], @(x)(isempty(x) || isstruct(x)));
    p.addParameter('fontSize', 16, @isscalar);

    p.parse(varargin{:});

    maxVisualizedRFs = p.Results.maxVisualizedRFs;
    onlyForRGCwithIndex = p.Results.onlyForRGCwithIndex;
    generateVideo = p.Results.generateVideo;
    videoFileName = p.Results.videoFileName;
    eccentricityCrossHairs = p.Results.withEccentricityCrossHairs;
    thePSFData = p.Results.withPSFData;
    fontSize = p.Results.fontSize;

    if (isempty(obj.rgcRFcenterConePoolingMatrix))
        fprintf(2,'The center and surround cone pooling matrices have not yet been set (no center/surround RF and no overlap) !!\n');
        return;
    end

    % Compute center of mosaic
    mRGCmosaicCenterDegs = mean(obj.rgcRFpositionsDegs,1);
   
    if (isempty(onlyForRGCwithIndex))
        % Sort RGCs according to their distance from the mosaic center
        ecc = sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2);
        [~,sortedRGCindices] = sort(ecc, 'ascend');
    
        if (maxVisualizedRFs < 0)
            % If negative half of the visualizedRFs will be from the mosaic
            % center and the other half from its 4 corners
            maxVisualizedRFs = -maxVisualizedRFs;
            m1 = min([floor(maxVisualizedRFs/2) numel(sortedRGCindices)]);
            sortedRGCindices = unique([sortedRGCindices(1:m1); sortedRGCindices(end:-1:numel(sortedRGCindices)-m1)]);
        end
    else
        sortedRGCindices = onlyForRGCwithIndex;
    end


    if (generateVideo)
        if (isempty(videoFileName))
            videoFileName = 'RetinalRFs.mp4';
        end
        videoOBJ = VideoWriter('RetinalRFs.mp4', 'MPEG-4');
        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();
    end

    hFig = figure(4000); clf;
    set(hFig, 'Position', [10 10 1450 900], 'Color', [1 1 1]);

    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
       'rowsNum', 2, ...
       'colsNum', 3, ...
       'heightMargin',  0.12, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.08, ...
       'topMargin',      0.02);

    for i = 1:size(subplotPosVectors,1)
        for j = 1:size(subplotPosVectors,1)
            allAxes{i,j} = subplot('Position', subplotPosVectors(i,j).v);
        end
    end

    for iSortedRGCindex = 1:numel(sortedRGCindices)

        if (iSortedRGCindex > maxVisualizedRFs)
            continue;
        end

        % Retrieve the RGCindex
        iRGC = sortedRGCindices(iSortedRGCindex);

        % Retrieve this cell's center cone indices and weights
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, iRGC)));
        indicesOfCenterCones = find(connectivityVector > 0.0001);
        weightsOfCenterCones = connectivityVector(indicesOfCenterCones);

        % Retrieve the correct RTVFTobj based on this cells position and
        % #of center cones. For now only checking the centerConesNum
%         iObj = find(...
%             (obj.theConesNumPooledByTheRFcenterGrid == numel(indicesOfCenterCones)) & ...  % match the conesNum in the center
%             (distancesToEccGrid == min(distancesToEccGrid)) ...                     % match to the closest eccGrid position
%             )
        iObj = find(...
            (obj.theConesNumPooledByTheRFcenterGrid == numel(indicesOfCenterCones)) ...  % match the conesNum in the center
            );
        theRTVFTobj = obj.theRetinaToVisualFieldTransformerOBJList{iObj};

        % Extract the model constants
        modelConstants = theRTVFTobj.rfComputeStruct.modelConstants;

        connectivityVector = full(obj.rgcRFsurroundConePoolingMatrix(:, iRGC));
        indicesOfSurroundCones = find(connectivityVector > 0.0001);
        weightsOfSurroundCones = connectivityVector(indicesOfSurroundCones);

        % Find the cones that contribute to 95% of the total surround
        includedSurroundPercentage = 0.9;
        totalSurroundWeight = sum(weightsOfSurroundCones(:));
        surroundConeInclusionWeight = includedSurroundPercentage*totalSurroundWeight;

        [~,idx] = sort(weightsOfSurroundCones, 'descend');

        theCumSumWeight = 0.0;
        id = 1; keepGoing = true;
        while (id < numel(idx)) &&( keepGoing)
            theCumSumWeight = theCumSumWeight + weightsOfSurroundCones(idx(id));
            if (theCumSumWeight >= surroundConeInclusionWeight)                
                keepGoing = false;
            else
                indicesOfIncludedSurroundCones = indicesOfSurroundCones(idx(1:id));
            end
            id = id + 1;
        end


        % Spatial support for retinal cone maps
        spatialSupportDegs(:,1) = -0.2:0.001:0.2;
        spatialSupportDegs(:,2) = spatialSupportDegs(:,1);

        % Compute theRetinalRFcenterConeMap
        theConePositionsDegs = obj.inputConeMosaic.coneRFpositionsDegs(indicesOfCenterCones,:);
        rfCenterPositionDegs = mean(theConePositionsDegs,1);
        theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
        theRetinalRFcenterConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
            obj.inputConeMosaic,...
            theConePositionsDegs, ...
            weightsOfCenterCones, ...
            spatialSupportDegs);

        % Compute theRetinalRFsurroundConeMap
        theConePositionsDegs = obj.inputConeMosaic.coneRFpositionsDegs(indicesOfSurroundCones,:);
        theConePositionsDegs = bsxfun(@minus, theConePositionsDegs, rfCenterPositionDegs);
        theRetinalRFsurroundConeMap = RetinaToVisualFieldTransformer.retinalSubregionConeMapFromPooledConeInputs(...
            obj.inputConeMosaic,...
            theConePositionsDegs, ...
            weightsOfSurroundCones, ...
            spatialSupportDegs);

        % Included surround cones
        theSurroundConePositions = bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsDegs(indicesOfIncludedSurroundCones,:), rfCenterPositionDegs);
        theSurroundConeDistancesFromTheRFcenter = bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsDegs(indicesOfIncludedSurroundCones,:), rfCenterPositionDegs);
        surroundRadiusDegs = max(sqrt(sum(theSurroundConeDistancesFromTheRFcenter.^2,2)));
        surroundOutline.x = rfCenterPositionDegs(1)+surroundRadiusDegs*cosd(0:10:360);
        surroundOutline.y = rfCenterPositionDegs(2)+surroundRadiusDegs*sind(0:10:360);

        % Spatial profiles along x, y
        xCenterProfile = sum(theRetinalRFcenterConeMap,1);
        xSurroundProfile = sum(theRetinalRFsurroundConeMap,1);
        yCenterProfile = sum(theRetinalRFcenterConeMap,2);
        ySurroundProfile = sum(theRetinalRFsurroundConeMap,2);

       
        xSupportDegs = spatialSupportDegs(:,1) + rfCenterPositionDegs(1);
        ySupportDegs = spatialSupportDegs(:,2) + rfCenterPositionDegs(2);

        profileRange = [-1/4 1] * 0.15; % [-max([max(xSurroundProfile(:)) max(xSurroundProfile(:))]) max([max(xCenterProfile(:)) max(yCenterProfile(:))])]
        profileTicks = profileRange(1): profileRange(2)/5: profileRange(2);

        mRGCmosaicXLims = [min(obj.rgcRFpositionsDegs(:,1)) max(obj.rgcRFpositionsDegs(:,1))];
        mRGCmosaicYLims = [min(obj.rgcRFpositionsDegs(:,2)) max(obj.rgcRFpositionsDegs(:,2))];

        spatialTicksX = linspace(mRGCmosaicXLims(1), mRGCmosaicXLims(end), 7);
        spatialTicksY = linspace(mRGCmosaicYLims(1), mRGCmosaicYLims(end), 7);

        if (eccentricityCrossHairs)
            mosaicCenterDegs = mean(obj.rgcRFpositionsDegs,1);
            radius = spatialTicksX(2)-spatialTicksX(1);
            for iRadius = 1:3
                crossHairsOutline(iRadius,1,:) = mosaicCenterDegs(1) + iRadius*radius*cosd(0:4:360);
                crossHairsOutline(iRadius,2,:) = mosaicCenterDegs(2) + iRadius*radius*sind(0:4:360);
            end

        end


        % The position of this RGC in the mosaic
        ax = subplot('Position', subplotPosVectors(1,1).v);
        hold(ax, 'off');
        plot(ax,obj.rgcRFpositionsDegs(iRGC,1), obj.rgcRFpositionsDegs(iRGC,2), 'kx', 'MarkerSize', 16, 'LineWidth', 1.5);
        axis(ax, 'equal');
        grid(ax, 'on');
        xlabel(ax, 'space, x (degs)')
        ylabel(ax, 'space, y (degs)')
        set(ax, 'XLim', mRGCmosaicXLims, 'YLim', mRGCmosaicYLims);
        set(ax, 'XTick', spatialTicksX, 'XTickLabel', sprintf('%2.2f\n', spatialTicksX));
        set(ax, 'YTick', spatialTicksY, 'YTickLabel', sprintf('%2.2f\n', spatialTicksY));
        set(ax, 'FontSize', fontSize);

        inputConeTypes = obj.inputConeMosaic.coneTypes(indicesOfCenterCones);
        if (numel(inputConeTypes) <6)
            if (numel(inputConeTypes) == 1)
                coneInfoString = sprintf('input cone: ');
            else
                coneInfoString = sprintf('input cones: ', numel(inputConeTypes));
            end
            for iCone = 1:numel(inputConeTypes)
                switch (inputConeTypes(iCone))
                    case cMosaic.LCONE_ID
                        coneInfoString = sprintf('%s L', coneInfoString);
                    case cMosaic.MCONE_ID
                        coneInfoString = sprintf('%s M', coneInfoString);
                    case cMosaic.SCONE_ID
                        coneInfoString = sprintf('%s S', coneInfoString);
                end
            end
        else
            coneInfoString = sprintf('%d input cones', numel(inputConeTypes));
        end
        coneInfoString = sprintf('%s (net weight: %2.2f)', coneInfoString, sum(weightsOfCenterCones));

        % The center cone map
        cMapEntries = 1024;
        cMap = (brewermap(cMapEntries, '*greys')).^0.5;
        ax = subplot('Position', subplotPosVectors(1,2).v);
        hold(ax, 'off');
        imagesc(ax,xSupportDegs, ySupportDegs, theRetinalRFcenterConeMap);
        hold(ax, 'on');

        if (isfield(obj.inputConeMosaic.coneApertureModifiers, 'shape') && (strcmp(obj.inputConeMosaic.coneApertureModifiers.shape, 'Gaussian')))
                gaussianSigma = obj.inputConeMosaic.coneApertureModifiers.sigma;
                visualizedApertures = 3*sqrt(2)*gaussianSigma * obj.inputConeMosaic.coneApertureDiametersDegs(indicesOfCenterCones);
            else
                fprintf(2,'cone aperture is not Gaussian, so cannot visualize characteristic radius. Visualizing the diameter\n');
                visualizedApertures = obj.inputConeMosaic.coneApertureDiametersDegs(indicesOfCenterCones);
        end

        for iCone = 1:numel(indicesOfCenterCones)
            xyCenter = obj.inputConeMosaic.coneRFpositionsDegs(indicesOfCenterCones(iCone),:);
            xConeOutline = xyCenter(1) + visualizedApertures(iCone)*0.5*cosd(0:20:360);
            yConeOutline = xyCenter(2) + visualizedApertures(iCone)*0.5*sind(0:20:360);
            switch (obj.inputConeMosaic.coneTypes(indicesOfCenterCones(iCone)))
                case cMosaic.LCONE_ID
                    coneColor = obj.inputConeMosaic.lConeColor;
                case cMosaic.MCONE_ID
                    coneColor = obj.inputConeMosaic.mConeColor*0.5;
                case cMosaic.SCONE_ID
                    coneColor = obj.inputConeMosaic.sConeColor;
                case cMosaic.KCONE_ID
                    coneColor = obj.inputConeMosaic.kConeColor;
            end
            plot(ax,xConeOutline,yConeOutline,'k-', 'Color', coneColor, 'LineWidth', 1.0);
        end

        if (~isempty(thePSFData))
            mosaicCenterDegs = obj.eccentricityDegs; %mean(obj.rgcRFpositionsDegs,1);
            xSupportDegsForPSF = thePSFData.psfSupportXdegs + mosaicCenterDegs(1);
            ySupportDegsForPSF = thePSFData.psfSupportYdegs + mosaicCenterDegs(2);
            alpha = 0.0;
            contourLineColor = [0 0 0.5];
            contourLineWidth = 0.5;
            cMosaic.semiTransparentContourPlot(ax, xSupportDegsForPSF, ySupportDegsForPSF, ...
                thePSFData.vLambdaWeightedData/max(thePSFData.vLambdaWeightedData(:))*max(abs(theRetinalRFcenterConeMap(:))), ...
                (-1.0:0.1:1.0)*max(abs(theRetinalRFcenterConeMap(:))), cMap, alpha, contourLineColor, ...
                'lineWidth', contourLineWidth, ...
                'edgeAlpha', 0.7);
            
        end

        if (eccentricityCrossHairs)
            for iRadius = 1:3
                plot(ax, squeeze(crossHairsOutline(iRadius,1,:)), ...
                      squeeze(crossHairsOutline(iRadius,2,:)), 'k-', 'LineWidth', 0.7, ...
                      'Color', [0.5 0.5 0.5]);
            end
        end

        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', mRGCmosaicXLims, 'YLim', mRGCmosaicYLims);
        colormap(ax, cMap);
        set(ax, 'Color', cMap(cMapEntries/2,:));
        set(ax, 'CLim', 0.7*max(theRetinalRFcenterConeMap(:))*[-1 1]);
        
        set(ax, 'XTick', spatialTicksX, 'XTickLabel', sprintf('%2.2f\n', spatialTicksX));
        set(ax, 'YTick', spatialTicksY, 'YTickLabel', sprintf('%2.2f\n', spatialTicksY));
        ylabel(ax, 'space, y (degs)');
        grid(ax, 'on');
        xtickangle(ax, 90);
        set(ax, 'FontSize', fontSize);
        title(ax,coneInfoString, 'FontWeight', 'Normal');
        

        % The center and surround y-profiles
        ax = subplot('Position', subplotPosVectors(1,3).v); cla(ax);
        hold(ax, 'off');
        shadedAreaBetweenTwoLines(ax, yCenterProfile', ySupportDegs', ...
            ySupportDegs'*0, cMap(512+256,:), 'none', 0.6, 1.5, '-');
        hold(ax, 'on');
        shadedAreaBetweenTwoLines(ax, -ySurroundProfile', ySupportDegs', ...
            ySupportDegs'*0, cMap(512-256,:), 'none', 0.6, 1.5, '-');

        plot(ax,yCenterProfile-ySurroundProfile, ySupportDegs, 'k-', 'LineWidth', 1.0);
        plot(ax,mRGCmosaicYLims*0, mRGCmosaicYLims, 'k--', 'LineWidth', 0.5);
        set(ax, 'YLim',  mRGCmosaicYLims);
        set(ax, 'XLim', profileRange);
        set(ax, 'XTick', profileTicks, 'XTickLabel', sprintf('%2.1f\n', profileTicks/max(profileTicks)));
        set(ax, 'YTick', spatialTicksY, 'YTickLabel', sprintf('%2.2f\n', spatialTicksY));
        grid(ax, 'on');
        axis(ax, 'square');
        set(ax, 'Xdir', 'reverse', 'YAxisLocation', 'right');
        set(ax, 'FontSize', fontSize);
        xlabel(ax, 'gain');
        xtickangle(ax, 0);
        title(ax,'retinal X- and Y-line weighting profiles', 'FontWeight', 'Normal');


        % The fitted model values
        ax = subplot('Position', subplotPosVectors(2,1).v);
        hold(ax, 'off');
        RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, theRTVFTobj.rfComputeStruct.retinalConePoolingParams);

        % The surround cone map
        ax = subplot('Position', subplotPosVectors(2,2).v);
        hold(ax, 'off');
        imagesc(ax,xSupportDegs, ySupportDegs, -theRetinalRFsurroundConeMap);
        hold(ax, 'on');

        if (isfield(obj.inputConeMosaic.coneApertureModifiers, 'shape') && (strcmp(obj.inputConeMosaic.coneApertureModifiers.shape, 'Gaussian')))
                gaussianSigma = obj.inputConeMosaic.coneApertureModifiers.sigma;
                visualizedApertures = 3*sqrt(2)*gaussianSigma * obj.inputConeMosaic.coneApertureDiametersDegs(indicesOfIncludedSurroundCones);
            else
                fprintf(2,'cone aperture is not Gaussian, so cannot visualize characteristic radius. Visualizing the diameter\n');
                visualizedApertures = obj.inputConeMosaic.coneApertureDiametersDegs(indicesOfIncludedSurroundCones);
        end

        for iCone = 1:numel(indicesOfIncludedSurroundCones)
            xyCenter = obj.inputConeMosaic.coneRFpositionsDegs(indicesOfIncludedSurroundCones(iCone),:);
            xConeOutline = xyCenter(1) + visualizedApertures(iCone)*0.5*cosd(0:20:360);
            yConeOutline = xyCenter(2) + visualizedApertures(iCone)*0.5*sind(0:20:360);
            switch (obj.inputConeMosaic.coneTypes(indicesOfIncludedSurroundCones(iCone)))
                case cMosaic.LCONE_ID
                    coneColor = obj.inputConeMosaic.lConeColor;
                case cMosaic.MCONE_ID
                    coneColor = obj.inputConeMosaic.mConeColor*0.5;
                case cMosaic.SCONE_ID
                    coneColor = obj.inputConeMosaic.sConeColor;
                case cMosaic.KCONE_ID
                    coneColor = obj.inputConeMosaic.kConeColor;
            end
            plot(ax,xConeOutline,yConeOutline,'k-', 'Color', coneColor, 'LineWidth', 1.0);
        end

        %plot(ax, surroundOutline.x, surroundOutline.y, 'k--', 'LineWidth', 1.0);
        if (eccentricityCrossHairs)
            for iRadius = 1:3
                plot(ax, squeeze(crossHairsOutline(iRadius,1,:)), ...
                      squeeze( crossHairsOutline(iRadius,2,:)), 'k-', 'LineWidth', 0.7, ...
                      'Color', [0.5 0.5 0.5]);
            end
        end
        
        hold(ax, 'off');
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', mRGCmosaicXLims, 'YLim', mRGCmosaicYLims)
        set(ax, 'CLim', 0.05*max(theRetinalRFsurroundConeMap(:))*[-1 1]);
        set(ax, 'Color', cMap(cMapEntries/2,:));
        xlabel(ax, 'space, x (degs)')
        ylabel(ax, 'space, y (degs)')
        set(ax, 'XTick', spatialTicksX, 'XTickLabel', sprintf('%2.2f\n', spatialTicksX));
        set(ax, 'YTick', spatialTicksY, 'YTickLabel', sprintf('%2.2f\n', spatialTicksY));
        grid(ax, 'on');
        xtickangle(ax, 90);
        set(ax, 'FontSize', fontSize);
        title(ax,sprintf('%d/%d surround cones accounting for\n%2.0f%% of the net surround weight (%2.3f)', ...
            numel(indicesOfIncludedSurroundCones), numel(indicesOfSurroundCones), includedSurroundPercentage*100, sum(weightsOfSurroundCones)), 'FontWeight', 'Normal');
        colormap(ax, cMap);

        % The center and surround x-profiles
        ax = subplot('Position', subplotPosVectors(2,3).v); cla(ax);
        hold(ax, 'off');
        shadedAreaBetweenTwoLines(ax, xSupportDegs', xCenterProfile, ...
            xCenterProfile*0, cMap(512+256,:), 'none', 0.6, 1.5, '-');
        hold(ax, 'on');
        shadedAreaBetweenTwoLines(ax, xSupportDegs', -xSurroundProfile, ...
            -xSurroundProfile*0, cMap(512-256,:), 'none', 0.6, 1.5, '-');
        plot(ax,xSupportDegs, xCenterProfile-xSurroundProfile, 'k-', 'LineWidth', 1.0);
        plot(ax,mRGCmosaicXLims, mRGCmosaicXLims*0, 'k--', 'LineWidth', 0.5);
        set(ax, 'XLim', mRGCmosaicXLims);
        set(ax, 'YLim', profileRange);
        set(ax, 'YTick', profileTicks, 'YTickLabel', sprintf('%2.1f\n', profileTicks/max(profileTicks)));
        set(ax, 'XTick', spatialTicksX, 'XTickLabel', sprintf('%2.2f\n', spatialTicksX));
        set(ax, 'FontSize', fontSize);
        xtickangle(ax, 90);
        grid(ax, 'on');
        xlabel(ax, 'space, x (degs)')
        ylabel(ax, 'gain');
        axis(ax, 'square');
        title(ax,sprintf('retinal S/C int. sens. ratio = %2.3f', sum(weightsOfSurroundCones)/sum(weightsOfCenterCones)), 'FontWeight', 'Normal');

        drawnow;
        if (generateVideo)
            videoOBJ.writeVideo(getframe(hFig));
        end
    end

    if (generateVideo)
        videoOBJ.close();
    end
end

function shadedAreaBetweenTwoLines(ax,x,y1, y2, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    x = [x  x(end)  fliplr(x)  x(1)];
    y = [y1 y2(end) fliplr(y2) y2(1)];
    px = reshape(x, [1 numel(x)]);
    py = reshape(y, [1 numel(y)]);
    pz = -10*eps*ones(size(py)); 
    patch(ax,px,py,pz,'FaceColor',faceColor,'EdgeColor', edgeColor, ...
        'FaceAlpha', faceAlpha, 'LineWidth', lineWidth, 'LineStyle', lineStyle);
end
