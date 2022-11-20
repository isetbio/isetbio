function [hFig, allAxes] = visualizeSpatialRFs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('maxVisualizedRFs', 18, @(x)(isempty(x) || isscalar(x)));
    p.addParameter('onlyForRGCwithIndex', [], @(x)(isempty(x) || isscalar(x)));
    p.addParameter('generateVideo', false, @islogical);

    p.parse(varargin{:});

    maxVisualizedRFs = p.Results.maxVisualizedRFs;
    onlyForRGCwithIndex = p.Results.onlyForRGCwithIndex;
    generateVideo = p.Results.generateVideo;

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
       'heightMargin',  0.1, ...
       'widthMargin',    0.05, ...
       'leftMargin',     0.03, ...
       'rightMargin',    0.01, ...
       'bottomMargin',   0.07, ...
       'topMargin',      0.03);

    for i = 1:size(subplotPosVectors,1)
        for j = 1:size(subplotPosVectors,1)
            allAxes{i,j} = subplot('Position', subplotPosVectors(i,j).v);
        end
    end
    
    % font size to use
    fontSize = 16;

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

        % Surround cones with weights > 0.1 x peak surround
        w = weightsOfSurroundCones;
        idx = find(w>=0.1*max(w));
        theSurroundConePositions = theConePositionsDegs(idx,:);
        surroundRadiusDegs = prctile(sqrt(sum(theSurroundConePositions.^2,2)), 95);
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
        ax = subplot('Position', subplotPosVectors(1,2).v);
        hold(ax, 'off');
        cMap = brewermap(1024, '*RdBu');
        imagesc(ax,xSupportDegs, ySupportDegs, theRetinalRFcenterConeMap);
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', mRGCmosaicXLims, 'YLim', mRGCmosaicYLims)
        set(ax, 'CLim', 0.7*max(theRetinalRFcenterConeMap(:))*[-1 1], 'Color', cMap(512,:));
        set(ax, 'XTick', spatialTicksX, 'XTickLabel', sprintf('%2.2f\n', spatialTicksX));
        set(ax, 'YTick', spatialTicksY, 'YTickLabel', sprintf('%2.2f\n', spatialTicksY));
        ylabel(ax, 'space, y (degs)');
        grid(ax, 'on');
        xtickangle(ax, 90);
        set(ax, 'FontSize', fontSize);
        

        title(ax,coneInfoString);
        colormap(ax, brewermap(1024, '*RdBu'));

        % The center and surround y-profiles
        ax = subplot('Position', subplotPosVectors(1,3).v); cla(ax);
        hold(ax, 'off');
        shadedAreaBetweenTwoLines(ax, yCenterProfile', ySupportDegs', ...
            ySupportDegs'*0, [1 0.5 0.5], 'none', 0.6, 1.5, '-');
        hold(ax, 'on');
        shadedAreaBetweenTwoLines(ax, -ySurroundProfile', ySupportDegs', ...
            ySupportDegs'*0, [0.5 0.5 1.0], 'none', 0.6, 1.5, '-');

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
        title(ax,sprintf('S/C = %2.3f', sum(weightsOfSurroundCones)/sum(weightsOfCenterCones)));


         % The fitted model values
        ax = subplot('Position', subplotPosVectors(2,1).v);
        hold(ax, 'off');
        RetinaToVisualFieldTransformer.visualizeRetinalSurroundModelParametersAndRanges(ax, theRTVFTobj.rfComputeStruct.retinalConePoolingParams);

        % The surround cone map
        ax = subplot('Position', subplotPosVectors(2,2).v);
        hold(ax, 'off');
        imagesc(ax,xSupportDegs, ySupportDegs, -theRetinalRFsurroundConeMap);
        hold(ax, 'on');
        plot(ax,rfCenterPositionDegs(1)+theSurroundConePositions(:,1), rfCenterPositionDegs(2)+theSurroundConePositions(:,2), 'ko');
        plot(ax, surroundOutline.x, surroundOutline.y, 'k--', 'LineWidth', 1.0);
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XLim', mRGCmosaicXLims, 'YLim', mRGCmosaicYLims)
        set(ax, 'CLim', 0.7*max(theRetinalRFsurroundConeMap(:))*[-1 1], 'Color', cMap(512,:));
        xlabel(ax, 'space, x (degs)')
        ylabel(ax, 'space, y (degs)')
        set(ax, 'XTick', spatialTicksX, 'XTickLabel', sprintf('%2.2f\n', spatialTicksX));
        set(ax, 'YTick', spatialTicksY, 'YTickLabel', sprintf('%2.2f\n', spatialTicksY));
        grid(ax, 'on');
        xtickangle(ax, 90);
        set(ax, 'FontSize', fontSize);
        title(ax,sprintf('%d surround cones (net weight: %2.2f)', numel(indicesOfSurroundCones), sum(weightsOfSurroundCones)));
        colormap(ax, brewermap(1024, '*RdBu'));

        % The center and surround x-profiles
        ax = subplot('Position', subplotPosVectors(2,3).v); cla(ax);
        hold(ax, 'off');
        shadedAreaBetweenTwoLines(ax, xSupportDegs', xCenterProfile, ...
            xCenterProfile*0, [1 0.5 0.5], 'none', 0.6, 1.5, '-');
        hold(ax, 'on');
        shadedAreaBetweenTwoLines(ax, xSupportDegs', -xSurroundProfile, ...
            -xSurroundProfile*0, [0.5 0.5 1], 'none', 0.6, 1.5, '-');
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
        title(sprintf('S/C = %2.3f', sum(weightsOfSurroundCones)/sum(weightsOfCenterCones)));

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
