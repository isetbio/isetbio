classdef osWindow < handle
    
    properties (Dependent)
        oi
        sensor
        os
    end
    
    properties (Access = private)
        
        oiPrivate;
        osPrivate;
        sensorPrivate; 
        
        % figure handle
        hFig;
        
        lastFigWidth = 0;
        lastFigHeight = 0;
        
        % struct containing all the figure axes
        axesStruct;
        
        % the time slider uicontrol
        timeSlider;
        
        % Optical image - related properties
        % - struct with handles to overlay plots for the optical image
        opticalImageOverlayPlots;
        
        % - RGB rendering of the optical image
        opticalImageRGBrendering;
        opticalImageRGBrenderingFullRes;
        
        % - X- and Y-axis for optical image (in microns)
        opticalImageXdata;
        opticalImageYdata;
        opticalImageXgrid;
        opticalImageYgrid;
        
        % Sensor - related properties
        sensorOutlineInMicrons;
        sensorSizeInMicrons
        sensorPositionsInMicrons;
        
        % - X- and Y-axis for sensor (in microns)
        % - struct with handles to overlay plots for the sensor view image
        sensorViewOverlayPlots;
        sensorViewXdata;
        sensorViewYdata;
        sensorViewOpticalImage;
        sensorXsamplingGrid;
        sensorYsamplingGrid;
        
        % Outer Segment - related properties
        % - struct with handles to overlay plots for the outersegment
        outerSegmentOverlayPlots;
         
        % - XYT and XT outersegment responses
        outerSegmentDisplayedCurrentRange;
        outerSegmentXYTCurrent;
        outerSegmentXTCurrent;
        outerSegmentResponseTimeData;
        outerSegmentResponseXYpositionData;
        outerSegmentResponseXpositionData;
        outerSegmentResponseYpositionData;
        outerSegmentResponseHiResXpositionData;
        outerSegmentResponseHiResYpositionData;
        outerSegmentResponseX;
        outerSegmentResponseY;
        outerSegmentResponseHiResX;
        outerSegmentResponseHiResY;
        outerSegmentSpatiallyInterpolated2DResponseMap;
        outerSegmentXYResponseInterpolatingKernel;
    end
    
    % Public API
    methods
        % Constructor
        function obj = osWindow(figureNo, figureTitle, os, sensor, oi)
            
            obj.init(figureNo, figureTitle); 
            
            obj.oi = oi;
            obj.sensor = sensor;
            obj.os = os;
            
            % whenever we set a new oi, sensor, we re-generate the different figure handles
            obj.generateAxesAndControls();
            obj.initOpticalImageDisplay();
            obj.initSensorViewDisplay();
            obj.initOuterSegmentResponseDisplays();
        end % constructor

        function set.os(obj, os)
            % generate our private copy of the outer segment
            obj.osPrivate = os;
            
            % get the full XYT response
            obj.outerSegmentXYTCurrent = osGet(obj.osPrivate, 'ConeCurrentSignal');
            
            % subtract baseline response at t = 0;
            obj.outerSegmentXYTCurrent = bsxfun(@minus, obj.outerSegmentXYTCurrent, squeeze(obj.outerSegmentXYTCurrent(:,:,1)));
            
            % XT response
            obj.outerSegmentXTCurrent  = reshape(obj.outerSegmentXYTCurrent, [size(obj.outerSegmentXYTCurrent,1)*size(obj.outerSegmentXYTCurrent,2), size(obj.outerSegmentXYTCurrent,3)]);
           
            % Displayed current range
            obj.outerSegmentDisplayedCurrentRange = [0 max(obj.outerSegmentXYTCurrent(:))];
                
            obj.outerSegmentResponseTimeData = (0:size(obj.outerSegmentXYTCurrent,3)-1)/size(obj.outerSegmentXYTCurrent,3)*sensorGet(obj.sensorPrivate, 'total time');
            obj.outerSegmentResponseXpositionData = 1:size(obj.outerSegmentXYTCurrent,2);
            obj.outerSegmentResponseYpositionData = 1:size(obj.outerSegmentXYTCurrent,1);
            delta = 0.10; sigma = 4.0*delta;
            obj.outerSegmentResponseHiResXpositionData = 0:delta:(obj.outerSegmentResponseXpositionData(end)+1);
            obj.outerSegmentResponseHiResYpositionData = 0:delta:(obj.outerSegmentResponseYpositionData(end)+1);
            [X,Y] = meshgrid(delta*(-10:10), delta*(-10:10)); 
            obj.outerSegmentXYResponseInterpolatingKernel = (exp(-(X/sigma).^2) .* exp(-(Y/sigma).^2)).^0.99;
            obj.outerSegmentResponseXYpositionData = 1:size(obj.outerSegmentXYTCurrent,1)*size(obj.outerSegmentXYTCurrent,2);
            
            % XY response helper properties
            [obj.outerSegmentResponseX, obj.outerSegmentResponseY] = meshgrid(obj.outerSegmentResponseXpositionData, obj.outerSegmentResponseYpositionData);
            [obj.outerSegmentResponseHiResX, obj.outerSegmentResponseHiResY] = meshgrid(obj.outerSegmentResponseHiResXpositionData, obj.outerSegmentResponseHiResYpositionData);
                             
        end
        
        function set.sensor(obj, sensor)
            % generate our private copy of the outer segment
            obj.sensorPrivate = sensor;
            
            % compute sensor positions in microns
            sensorSampleSeparationInMicrons = sensorGet(obj.sensorPrivate,'pixel size','um');
            pos = sensorGet(obj.sensorPrivate,'positions');
            obj.sensorPositionsInMicrons(:,1) = -pos(:,1)*sensorSampleSeparationInMicrons(1);
            obj.sensorPositionsInMicrons(:,2) =  pos(:,2)*sensorSampleSeparationInMicrons(2);
            %obj.sensorPositionsInMicrons = -bsxfun(@times, sensorGet(obj.sensorPrivate,'positions'), sensorSampleSeparationInMicrons);
            
            % compute sensor cone sampling grid
            sensorRowsCols = sensorGet(obj.sensorPrivate, 'size');
            dx = sensorRowsCols(2) * sensorSampleSeparationInMicrons(2);
            dy = sensorRowsCols(1) * sensorSampleSeparationInMicrons(1);
            obj.sensorSizeInMicrons = [dx dy]; 
           
            [R,C] = meshgrid(1:sensorRowsCols(1), 1:sensorRowsCols(2));
            obj.sensorXsamplingGrid = (C(:)-0.5) * sensorSampleSeparationInMicrons(1);
            obj.sensorYsamplingGrid = (R(:)-0.5) * sensorSampleSeparationInMicrons(2);
            obj.sensorOutlineInMicrons(:,1) = [-1 -1 1 1 -1] * dx/2;
            obj.sensorOutlineInMicrons(:,2) = [-1 1 1 -1 -1] * dy/2;
        end
        
        function set.oi(obj, oi)
            % generate our private copy of the optical image
            obj.oiPrivate = oi;
            
            % generate RGB rendering of the optical image
            obj.opticalImageRGBrenderingFullRes = oiGet(obj.oiPrivate, 'rgb image');
            
            oiSpatialSupport = oiGet(obj.oiPrivate,'spatial support','microns');
            obj.opticalImageXgrid = squeeze(oiSpatialSupport(:,:,1)); 
            obj.opticalImageYgrid = squeeze(oiSpatialSupport(:,:,2));
            obj.opticalImageXdata = squeeze(obj.opticalImageXgrid(1,:));  % x-positions from 1st row
            obj.opticalImageYdata = squeeze(obj.opticalImageYgrid(:,1));  % y-positions from 1st col
    
            % subsample the optical image by a factor of 2 to speed up display
            k = 2;
            obj.opticalImageXdata = obj.opticalImageXdata(1):k:obj.opticalImageXdata(end);
            obj.opticalImageYdata = obj.opticalImageYdata(1):k:obj.opticalImageYdata(end);
            obj.opticalImageRGBrendering = obj.opticalImageRGBrenderingFullRes(1:k:end, 1:k:end,:);
        end
    
    end
        
    methods (Access = private)  
        
        function init(obj, figureNo, figureTitle)
            obj.hFig = figure(figureNo);
            aspectRatio = 800/1000;
            screenSize = get(0,'ScreenSize');
            screenSize(4) = screenSize(4)*0.85;
            set(obj.hFig, 'Name', figureTitle, 'Menubar', 'none', 'Toolbar', 'none', 'Color', [0.1 0.1 0.1], ...
                'Position',[10+rand(1,1)*300 10+rand(1,1)*100 screenSize(4)*aspectRatio screenSize(4)], ...
                'SizeChangedFcn', {@resizeOSwindow, obj, aspectRatio});
        end
        
        function initOuterSegmentResponseDisplays(obj)
            positionIndex = 10;
            currentTime = obj.outerSegmentResponseTimeData(positionIndex);
            N = numel(obj.outerSegmentResponseXYpositionData);
            
            % The XT plot
            cla(obj.axesStruct.outerSegmentXTresponseAxes);
            obj.outerSegmentOverlayPlots.p1 = imagesc(...
                'XData', obj.outerSegmentResponseTimeData, 'YData', obj.outerSegmentResponseXYpositionData, ...
                'CData', obj.outerSegmentXTCurrent, 'parent', obj.axesStruct.outerSegmentXTresponseAxes);
            % plot current time vertical line
            hold(obj.axesStruct.outerSegmentXTresponseAxes, 'on');
            obj.outerSegmentOverlayPlots.p2  = plot(obj.axesStruct.outerSegmentXTresponseAxes, currentTime*[1 1], [1 N] , 'y--', 'LineWidth', 2);
            hold(obj.axesStruct.outerSegmentXTresponseAxes, 'off');
            set(obj.axesStruct.outerSegmentXTresponseAxes, 'XTickLabel', {});
            axis(obj.axesStruct.outerSegmentXTresponseAxes, 'ij');
            c = colorbar(obj.axesStruct.outerSegmentXTresponseAxes, 'northoutside');
            set(c, 'Color', [0.5 0.5 0.5], 'FontSize', 9);
            colormap(obj.axesStruct.outerSegmentXTresponseAxes, 'bone');
            
            set(obj.axesStruct.outerSegmentXTresponseAxes, ...
                 'XLim', [obj.outerSegmentResponseTimeData(1) obj.outerSegmentResponseTimeData(end)], ...
                 'YLim', [obj.outerSegmentResponseXYpositionData(1) obj.outerSegmentResponseXYpositionData(end)], ...
                 'CLim', obj.outerSegmentDisplayedCurrentRange, ...
                 'XColor', [1 1 1], 'YColor', [1 1 1]);
            ylabel(obj.axesStruct.outerSegmentXTresponseAxes,  'cone #', 'FontSize', 12);
            set(obj.axesStruct.outerSegmentXTresponseAxes, 'FontSize', 12);

            % The traces plot
            cla(obj.axesStruct.outerSegmentTracesAxes);
            hold(obj.axesStruct.outerSegmentTracesAxes, 'on');
            
            [lConeIndex, mConeIndex, sConeIndex] = obj.findStrongestResponsingCone();
            
            coneTypes = sensorGet(obj.sensorPrivate, 'cone type');
            [rowBestLcone, colBestLcone] = ind2sub(size(coneTypes), lConeIndex);
            [rowBestMcone, colBestMcone] = ind2sub(size(coneTypes), mConeIndex);
            [rowBestScone, colBestScone] = ind2sub(size(coneTypes), sConeIndex);
            
            obj.outerSegmentOverlayPlots.p4 = plot(obj.axesStruct.outerSegmentTracesAxes, obj.outerSegmentResponseTimeData, obj.outerSegmentXTCurrent(lConeIndex,:), 'r-', 'Color', [1 0.2 0.4]);
            obj.outerSegmentOverlayPlots.p5 = plot(obj.axesStruct.outerSegmentTracesAxes, obj.outerSegmentResponseTimeData, obj.outerSegmentXTCurrent(mConeIndex,:), 'g-', 'Color', [0.2 1.0 0.5]);
            obj.outerSegmentOverlayPlots.p6 = plot(obj.axesStruct.outerSegmentTracesAxes, obj.outerSegmentResponseTimeData, obj.outerSegmentXTCurrent(sConeIndex,:), 'b-', 'Color', [0.8 0.3 1.0]);
            obj.outerSegmentOverlayPlots.p7 = plot(obj.axesStruct.outerSegmentTracesAxes, currentTime*[1 1], obj.outerSegmentDisplayedCurrentRange , 'y--', 'LineWidth', 2);
            hold(obj.axesStruct.outerSegmentTracesAxes, 'off');
            box(obj.axesStruct.outerSegmentTracesAxes, 'on');
            set(obj.axesStruct.outerSegmentTracesAxes, ...
                 'XLim', [obj.outerSegmentResponseTimeData(1) obj.outerSegmentResponseTimeData(end)], ...
                 'YLim', obj.outerSegmentDisplayedCurrentRange, ...
                 'XColor', [1 1 1], 'YColor', [1 1 1], 'Color', [0 0 0] ...
                 );
            xlabel(obj.axesStruct.outerSegmentTracesAxes, 'time (seconds)', 'FontSize', 12);
            ylabel(obj.axesStruct.outerSegmentTracesAxes, 'current (uAmps)', 'FontSize', 12);
            set(obj.axesStruct.outerSegmentTracesAxes, 'FontSize', 12);
            
            % The XY plot
            obj.computeSpatiallyInterpolatedOuterSegment2DResponseMap(positionIndex);

            cla(obj.axesStruct.outerSegmentXYresponseAxes);
            obj.outerSegmentOverlayPlots.p3 = ...
                 imagesc('XData', obj.outerSegmentResponseHiResXpositionData, ...
                         'YData', obj.outerSegmentResponseHiResYpositionData, ..., 
                         'CData', obj.outerSegmentSpatiallyInterpolated2DResponseMap, ...
                         'parent', obj.axesStruct.outerSegmentXYresponseAxes);
            hold(obj.axesStruct.outerSegmentXYresponseAxes, 'on');
            plot(obj.axesStruct.outerSegmentXYresponseAxes, obj.outerSegmentResponseXpositionData(colBestLcone), obj.outerSegmentResponseYpositionData(rowBestLcone), 'ro');
            plot(obj.axesStruct.outerSegmentXYresponseAxes, obj.outerSegmentResponseXpositionData(colBestMcone), obj.outerSegmentResponseYpositionData(rowBestMcone), 'go');
            plot(obj.axesStruct.outerSegmentXYresponseAxes, obj.outerSegmentResponseXpositionData(colBestScone), obj.outerSegmentResponseYpositionData(rowBestScone), 'bo');
            hold(obj.axesStruct.outerSegmentXYresponseAxes, 'off');
            set(obj.axesStruct.outerSegmentXYresponseAxes, ...
                'XLim', [0 max(obj.outerSegmentResponseHiResXpositionData)], ...
                'YLim', [0 max(obj.outerSegmentResponseHiResYpositionData)], ...
                'CLim', obj.outerSegmentDisplayedCurrentRange, ...
                'XColor', [1 1 1], 'YColor', [1 1 1]);
            
            axis(obj.axesStruct.outerSegmentXYresponseAxes,'ij'); axis(obj.axesStruct.outerSegmentXYresponseAxes,'equal');
            axis(obj.axesStruct.outerSegmentXYresponseAxes, 'off');
            set(obj.axesStruct.outerSegmentXYresponseAxes, 'FontSize', 12, 'XTick', [], 'YTick', []);
            colormap(obj.axesStruct.outerSegmentXYresponseAxes, 'bone');
            title(obj.axesStruct.outerSegmentXYresponseAxes, sprintf('t = %2.3f sec', obj.outerSegmentResponseTimeData(positionIndex)), 'Color', [0.9 0.7 0.1], 'FontSize', 12);
        end
        
        
        function updateOuterSegmentResponseDisplays(obj, kPos)    
            % update the time line in the XT plot
            currentTime = obj.outerSegmentResponseTimeData(kPos);
            set(obj.outerSegmentOverlayPlots.p2, 'XData', currentTime*[1 1]);
            
            % update the traces plot
            set(obj.outerSegmentOverlayPlots.p7, 'XData', currentTime*[1 1]);
            
            % update the XY plot
            obj.computeSpatiallyInterpolatedOuterSegment2DResponseMap(kPos);
            set(obj.outerSegmentOverlayPlots.p3, 'CData', obj.outerSegmentSpatiallyInterpolated2DResponseMap); 
            title(obj.axesStruct.outerSegmentXYresponseAxes, sprintf('t = %2.3f sec', obj.outerSegmentResponseTimeData(kPos)), 'Color', [0.9 0.7 0.1], 'FontSize', 12);
        end
        
        
        function initSensorViewDisplay(obj)
            positionIndex = 10;
            currentSensorPosition = squeeze(obj.sensorPositionsInMicrons(positionIndex,:));
            obj.findImagePixelsUnderSensor(currentSensorPosition);
            
            cla(obj.axesStruct.sensorViewAxes);
            obj.sensorViewOverlayPlots.p1 = image('XData', obj.sensorViewXdata, 'YData', obj.sensorViewYdata, 'CData', obj.sensorViewOpticalImage, 'parent', obj.axesStruct.sensorViewAxes);
            hold(obj.axesStruct.sensorViewAxes, 'on');
            xpos = currentSensorPosition(1) -obj.sensorSizeInMicrons(1)/2 +  obj.sensorXsamplingGrid;
            ypos = currentSensorPosition(2) -obj.sensorSizeInMicrons(2)/2 +  obj.sensorYsamplingGrid;
            obj.sensorViewOverlayPlots.p2  = plot(obj.axesStruct.sensorViewAxes, xpos, ypos, 'k.', 'MarkerSize', 12);
            hold(obj.axesStruct.sensorViewAxes, 'off');
            axis(obj.axesStruct.sensorViewAxes,'ij'); axis(obj.axesStruct.sensorViewAxes,'equal');
            
            set(obj.axesStruct.sensorViewAxes, ...
                 'XLim', round(currentSensorPosition(1) + obj.sensorSizeInMicrons(1)*0.55*[-1 1]), ...
                 'YLim', round(currentSensorPosition(2) + obj.sensorSizeInMicrons(2)*0.55*[-1 1]), ...
                 'XColor', [01 1 1], 'YColor', [1 1 1]);
            set(obj.axesStruct.sensorViewAxes, 'FontSize', 12);
            tickPositions = -2000:20:2000;
            set(obj.axesStruct.sensorViewAxes, 'XTick', tickPositions, 'YTick', tickPositions);
            box(obj.axesStruct.sensorViewAxes, 'on');
            axis(obj.axesStruct.sensorViewAxes, 'off');
        end
        
        function updateSensorViewDisplay(obj, kPos)
            currentSensorPosition = squeeze(obj.sensorPositionsInMicrons(kPos,:));
            obj.findImagePixelsUnderSensor(currentSensorPosition);
            set(obj.sensorViewOverlayPlots.p1, 'XData', obj.sensorViewXdata, 'YData', obj.sensorViewYdata, 'CData', obj.sensorViewOpticalImage);  
            xpos = currentSensorPosition(1) -obj.sensorSizeInMicrons(1)/2 +  obj.sensorXsamplingGrid;
            ypos = currentSensorPosition(2) -obj.sensorSizeInMicrons(2)/2 +  obj.sensorYsamplingGrid;
            set(obj.sensorViewOverlayPlots.p2, 'XData', xpos, 'YData', ypos); 
            set(obj.axesStruct.sensorViewAxes, ...
                 'XLim', round(currentSensorPosition(1) + obj.sensorSizeInMicrons(1)*0.55*[-1 1]), ...
                 'YLim', round(currentSensorPosition(2) + obj.sensorSizeInMicrons(2)*0.55*[-1 1]));
        end
        
        function findImagePixelsUnderSensor(obj,currentSensorPosition)
            % find image pixels falling within the sensor outline
            pixelIndices = find(...
                (obj.opticalImageXgrid >= currentSensorPosition(1) - obj.sensorSizeInMicrons(1)*0.6) & ...
                (obj.opticalImageXgrid <= currentSensorPosition(1) + obj.sensorSizeInMicrons(1)*0.6) & ...
                (obj.opticalImageYgrid >= currentSensorPosition(2) - obj.sensorSizeInMicrons(2)*0.6) & ...
                (obj.opticalImageYgrid <= currentSensorPosition(2) + obj.sensorSizeInMicrons(2)*0.6) );
            [rows, cols] = ind2sub(size(obj.opticalImageXgrid), pixelIndices);
            
            rowRange = min(rows):1:max(rows);
            colRange = min(cols):1:max(cols);
            obj.sensorViewOpticalImage = obj.opticalImageRGBrenderingFullRes(rowRange,colRange,:);
            xGridSubset = obj.opticalImageXgrid(rowRange, colRange);
            yGridSubset = obj.opticalImageYgrid(rowRange, colRange);
            obj.sensorViewXdata = squeeze(xGridSubset(1,:));
            obj.sensorViewYdata = squeeze(yGridSubset(:,1));
        end
        
        function [lConeIndex, mConeIndex, sConeIndex] = findStrongestResponsingCone(obj)
            coneTypes = sensorGet(obj.sensorPrivate, 'cone type');
            lConeIndices = find(coneTypes == 2);
            mConeIndices = find(coneTypes == 3);
            sConeIndices = find(coneTypes == 4);
            
            A = obj.outerSegmentXTCurrent(lConeIndices,:);
            [m,index] = max(A(:)); [row, col] = ind2sub(size(A),index);
            lConeIndex = lConeIndices(row);
            
            A = obj.outerSegmentXTCurrent(mConeIndices,:);
            [m,index] = max(A(:)); [row, col] = ind2sub(size(A),index);
            mConeIndex = mConeIndices(row);
            
            A = obj.outerSegmentXTCurrent(sConeIndices,:);
            [m,index] = max(A(:)); [row, col] = ind2sub(size(A),index);
            sConeIndex = sConeIndices(row);
        end
        
        function computeSpatiallyInterpolatedOuterSegment2DResponseMap(obj, kPos)  
            obj.outerSegmentSpatiallyInterpolated2DResponseMap = zeros(numel(obj.outerSegmentResponseHiResYpositionData), numel(obj.outerSegmentResponseHiResXpositionData));
            delta = obj.outerSegmentResponseHiResXpositionData(2)-obj.outerSegmentResponseHiResXpositionData(1); step = round(1.0/delta);
            obj.outerSegmentSpatiallyInterpolated2DResponseMap(1+(1:size(obj.outerSegmentXYTCurrent,1))*step, 1+(1:size(obj.outerSegmentXYTCurrent,2))*step) = squeeze(obj.outerSegmentXYTCurrent(:,:,kPos));
            obj.outerSegmentSpatiallyInterpolated2DResponseMap = conv2(obj.outerSegmentSpatiallyInterpolated2DResponseMap, obj.outerSegmentXYResponseInterpolatingKernel, 'same');
        end
        
        function initOpticalImageDisplay(obj)  
            positionIndex = 10;
            currentSensorPosition = squeeze(obj.sensorPositionsInMicrons(positionIndex,:));
            sensorPositionHistory = obj.sensorPositionsInMicrons(1:positionIndex,:);
            
            cla(obj.axesStruct.opticalImageAxes);
            image('XData', obj.opticalImageXdata, 'YData', obj.opticalImageYdata, 'CData', ...
                  obj.opticalImageRGBrendering, 'parent', obj.axesStruct.opticalImageAxes);
            hold(obj.axesStruct.opticalImageAxes, 'on');
            obj.opticalImageOverlayPlots.p1 = plot(obj.axesStruct.opticalImageAxes, currentSensorPosition(1) + obj.sensorOutlineInMicrons(:,1), currentSensorPosition(2) + obj.sensorOutlineInMicrons(:,2), 'r-', 'LineWidth', 2);
            obj.opticalImageOverlayPlots.p2 = plot(obj.axesStruct.opticalImageAxes, currentSensorPosition(1) + obj.sensorOutlineInMicrons(:,1), currentSensorPosition(2) + obj.sensorOutlineInMicrons(:,2), 'w-', 'LineWidth', 1);
            obj.opticalImageOverlayPlots.p3 = plot(obj.axesStruct.opticalImageAxes, sensorPositionHistory(:,1), sensorPositionHistory(:,2), 'k-', 'LineWidth', 1);
            hold(obj.axesStruct.opticalImageAxes, 'off');
            axis(obj.axesStruct.opticalImageAxes,'ij'); axis(obj.axesStruct.opticalImageAxes,'equal');
            set(obj.axesStruct.opticalImageAxes, 'XLim', 0.8*max(abs(obj.opticalImageXdata(:)))*[-1 1], 'YLim', .8*max(abs(obj.opticalImageYdata(:)))*[-1 1])
            set(obj.axesStruct.opticalImageAxes, 'XTick', [], 'YTick', []);
        end
        
        function updateOpticalImageDisplay(obj, kPos)
            currentSensorPosition = squeeze(obj.sensorPositionsInMicrons(kPos,:));
            sensorPositionHistory = obj.sensorPositionsInMicrons(1:kPos,:);
            set(obj.opticalImageOverlayPlots.p1, 'XData', currentSensorPosition(1) + obj.sensorOutlineInMicrons(:,1), 'YData', currentSensorPosition(2) + obj.sensorOutlineInMicrons(:,2));
            set(obj.opticalImageOverlayPlots.p2, 'XData', currentSensorPosition(1) + obj.sensorOutlineInMicrons(:,1), 'YData', currentSensorPosition(2) + obj.sensorOutlineInMicrons(:,2));
            set(obj.opticalImageOverlayPlots.p3, 'XData', sensorPositionHistory(:,1), 'YData', sensorPositionHistory(:,2));
        end
        
        function generateAxesAndControls(obj)  
            w = 800;
            h = 1000;
            imageWidthToHeightRatio = size(obj.opticalImageRGBrendering,2) / size(obj.opticalImageRGBrendering,1);
            
            leftMargin = 5/w;
            opticalImageWidth  = (w-10)/w;
            opticalImageHeight = (w-10)/imageWidthToHeightRatio/h;
            bottomMargin = (h-10)/h - opticalImageHeight + 5/h;
    
            % generate plot axes
            sensorViewWidth = 200/w; sensorViewHeight = 200/h; 
            spatiotemporalViewWidth = 500/w; spatiotemporalViewHeight = 200/h;
            spatialViewWidth  = 200/w; spatialViewHeight = 200/h;
            obj.axesStruct.opticalImageAxes = axes('parent',obj.hFig,'unit','normalized','position',[leftMargin bottomMargin opticalImageWidth opticalImageHeight], 'Color', [0 0 0]);
            obj.axesStruct.sensorViewAxes   = axes('parent',obj.hFig,'unit','normalized','position',[leftMargin+20/w bottomMargin+20/h sensorViewWidth sensorViewHeight], 'Color', [0 0 0]);
            
            % generate response axes
            obj.axesStruct.outerSegmentXTresponseAxes = axes('parent',obj.hFig, 'unit','normalized','position',[9*leftMargin bottomMargin-spatiotemporalViewHeight+5/h     spatiotemporalViewWidth spatiotemporalViewHeight], 'Color', [0 0 0]);
            obj.axesStruct.outerSegmentTracesAxes = axes('parent',    obj.hFig, 'unit','normalized','position',[9*leftMargin bottomMargin-1.5*spatiotemporalViewHeight-25/h spatiotemporalViewWidth spatiotemporalViewHeight/2], 'Color', [0 0 0]);
            
            % generate 2D instantaneous response axes
            positionVector = [5*leftMargin+50/w+spatiotemporalViewWidth bottomMargin-1.5*spatiotemporalViewHeight - 10/h spatialViewWidth spatialViewHeight];
            obj.axesStruct.outerSegmentXYresponseAxes = axes('parent',obj.hFig,'unit','normalized','position', positionVector, 'Color', [0 0 0]);
            
            % generate time slider
            timeSliderLeftMargin = leftMargin;
            timeSliderBottom = (5)/h;
            
            obj.timeSlider = uicontrol(...
                'Parent', obj.hFig,...
                'Style', 'slider',...
                'BackgroundColor', [0.8 0.7 0.2], ...
                'Min', 1, 'Max', size(obj.sensorPositionsInMicrons,1), 'Value', 1,...
                'Units', 'normalized',...
                'Position', [timeSliderLeftMargin, timeSliderBottom 0.99 0.012]);    
           
            % set the slider step
            set(obj.timeSlider, 'SliderStep', 1.0/((obj.timeSlider.Max-obj.timeSlider.Min)*10)*[1 1]);
            
            % set the callback
            addlistener(obj.timeSlider,'ContinuousValueChange', ...
                                      @(hFigure,eventdata) timeSliderCallback(obj.timeSlider,eventdata, obj));                          
        end
        
    end
end

% Callback for figure resizing
function resizeOSwindow(hObject,Event, obj, aspectRatio)
    posVector = get(hObject,'Position');
    width = posVector(3);
    height = width/aspectRatio;
    %Set(obj.hFig,'Position',[posVector(1:2) width height]);
end


% Callback for time slider
function timeSliderCallback(hObject,eventdata, obj)
    currentTimeBin = round(get(hObject,'Value'));
    obj.updateOpticalImageDisplay(currentTimeBin);
    obj.updateSensorViewDisplay(currentTimeBin);
    obj.updateOuterSegmentResponseDisplays(currentTimeBin);
end