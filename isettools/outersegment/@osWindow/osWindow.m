% Class for dynamic visualization of outer segment currents in response to
% eye movements scanning hyperspectral images.
% 
% For example usage see: t_osHyperSpectralSceneEyeScan
%
% 1/2016        NPC     Created.
% 2/24/2016     NPC     Added method for video generation, window resizing

classdef osWindow < handle
    
    properties (Dependent)
        scene
        oi
        sensor
        os 
    end
    
    properties (Access = private) 
        scenePrivate;
        oiPrivate;
        osPrivate;
        sensorPrivate; 
        zoomedInView;
        
        % figure handle
        hFig;
        
        % figure layout. Either 'horizontalLayout' or 'verticalLayout'
        figOrientation;
        
        % figure width to height aspect ratio
        aspectRatio;
        
        % struct containing all the figure axes
        axesStruct;
        
        % the sensor view uicontrols
        sensorViewText
        sensorViewSlider
        
        % the time slider uicontrol
        timeSlider;
        
        % the min, max displayed response
        minDisplayedReponseSlider
        maxDisplayedReponseSlider
        
        % Scene - related properties
        % - RGB rendering of the scene
        sceneRGBrenderingFullRes;
        
        % - X- and Y-axis for scene (in microns)
        sceneXdata;
        sceneYdata;
        sceneXgrid;
        sceneYgrid;
        
        
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
        sensorView;
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
        function obj = osWindow(figureNo, figureTitle, figOrientation, os, sensor, oi, scene)
            
            obj.init(figureNo, figureTitle, figOrientation); 
            
            obj.zoomedInView = 'optical image';
            
            obj.oi = oi;
            obj.scene = scene;
            obj.sensor = sensor;
            obj.os = os;
            
            % whenever we set a new oi, sensor, we re-generate the different figure handles
            obj.generateAxesAndControls(figOrientation);
            obj.initOpticalImageDisplay();
            obj.initSensorViewDisplay();
            obj.initSensorStimulusEncodingDisplay();
            obj.initOuterSegmentResponseDisplays();
        end % constructor

        function set.os(obj, os)
            % generate our private copy of the outer segment
            obj.osPrivate = os;
            
            % get the full XYT response
            obj.outerSegmentXYTCurrent = osGet(obj.osPrivate, 'ConeCurrentSignal');
            
            % subtract baseline response at t = 0;
            %obj.outerSegmentXYTCurrent = bsxfun(@minus, obj.outerSegmentXYTCurrent, squeeze(obj.outerSegmentXYTCurrent(:,:,1)));
            
            % XT response
            obj.outerSegmentXTCurrent  = reshape(obj.outerSegmentXYTCurrent, [size(obj.outerSegmentXYTCurrent,1)*size(obj.outerSegmentXYTCurrent,2), size(obj.outerSegmentXYTCurrent,3)]);
           
            % Displayed current range
            if isinf(max(obj.outerSegmentXYTCurrent(:)))
                errordlg('Overflow in computed OSX current. Use smaller time step', 'Fatal Error');
            end
            
            
            obj.outerSegmentResponseTimeData = (0:size(obj.outerSegmentXYTCurrent,3)-1)/size(obj.outerSegmentXYTCurrent,3)*sensorGet(obj.sensorPrivate, 'total time');
            obj.outerSegmentResponseXpositionData = 1:size(obj.outerSegmentXYTCurrent,2);
            obj.outerSegmentResponseYpositionData = 1:size(obj.outerSegmentXYTCurrent,1);
            delta = 0.05; sigma = 0.21;
            obj.outerSegmentResponseHiResXpositionData = 0:delta:(obj.outerSegmentResponseXpositionData(end)+1);
            obj.outerSegmentResponseHiResYpositionData = 0:delta:(obj.outerSegmentResponseYpositionData(end)+1);
            [X,Y] = meshgrid(delta*(-10:10), delta*(-10:10)); 
            obj.outerSegmentXYResponseInterpolatingKernel = (exp(-(X/sigma).^2) .* exp(-(Y/sigma).^2));
            obj.outerSegmentXYResponseInterpolatingKernel = obj.outerSegmentXYResponseInterpolatingKernel / max(obj.outerSegmentXYResponseInterpolatingKernel(:));
            obj.outerSegmentXYResponseInterpolatingKernel(obj.outerSegmentXYResponseInterpolatingKernel < 0.01) = 0;
            obj.outerSegmentXYResponseInterpolatingKernel = obj.outerSegmentXYResponseInterpolatingKernel.^0.5;
            
            
            obj.outerSegmentResponseXYpositionData = 1:size(obj.outerSegmentXYTCurrent,1)*size(obj.outerSegmentXYTCurrent,2);
            
            % XY response helper properties
            [obj.outerSegmentResponseX, obj.outerSegmentResponseY] = meshgrid(obj.outerSegmentResponseXpositionData, obj.outerSegmentResponseYpositionData);
            [obj.outerSegmentResponseHiResX, obj.outerSegmentResponseHiResY] = meshgrid(obj.outerSegmentResponseHiResXpositionData, obj.outerSegmentResponseHiResYpositionData);                 
        
            t = find(obj.outerSegmentResponseTimeData > 0.5);
            minTimIndex = t(1); 
            obj.outerSegmentDisplayedCurrentRange = [min(min(min(obj.outerSegmentXYTCurrent(:,:,minTimIndex:end)))) max(obj.outerSegmentXYTCurrent(:))];
            set(obj.minDisplayedReponseSlider, 'Value', obj.outerSegmentDisplayedCurrentRange(1));
            set(obj.maxDisplayedReponseSlider, 'Value', obj.outerSegmentDisplayedCurrentRange(2));
            
        end
        
        function set.sensor(obj, sensor)
            % generate our private copy of the outer segment
            obj.sensorPrivate = sensor;
            
            % compute sensor positions in microns
            sensorSampleSeparationInMicrons = sensorGet(obj.sensorPrivate,'pixel size','um');
            pos = sensorGet(obj.sensorPrivate,'positions');
            obj.sensorPositionsInMicrons(:,1) = -pos(:,1)*sensorSampleSeparationInMicrons(1);
            obj.sensorPositionsInMicrons(:,2) =  pos(:,2)*sensorSampleSeparationInMicrons(2);
            
            % compute sensor cone sampling grid
            sensorRowsCols = sensorGet(obj.sensorPrivate, 'size');
            dx = sensorRowsCols(2) * sensorSampleSeparationInMicrons(2);
            dy = sensorRowsCols(1) * sensorSampleSeparationInMicrons(1);
            obj.sensorSizeInMicrons = [dx dy];
           
            [R,C] = meshgrid(1:sensorRowsCols(1), 1:sensorRowsCols(2));
            R = R'; C = C';
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
    
            % subsample the displayed optical image by a factor of 2 to speed up display
            k = 2;
            obj.opticalImageXdata = obj.opticalImageXdata(1):k:obj.opticalImageXdata(end);
            obj.opticalImageYdata = obj.opticalImageYdata(1):k:obj.opticalImageYdata(end);
            obj.opticalImageRGBrendering = obj.opticalImageRGBrenderingFullRes(1:k:end, 1:k:end,:);
        end
        
        function set.scene(obj, scene)
            % generate our private copy of the scene
            obj.scenePrivate = scene;
            
            % generate RGB rendering of the optical image
            obj.sceneRGBrenderingFullRes = sceneGet(obj.scenePrivate, 'rgb image');
            
            % compute the retinal microns per degree conversion factor
            retinalMicronsPerPixel = oiGet(obj.oiPrivate,'wres','microns');
            degreesPerPixel = oiGet(obj.oiPrivate, 'angularresolution');
            retinalMicronsPerDegreeX = retinalMicronsPerPixel / degreesPerPixel(1);
            retinalMicronsPerDegreeY = retinalMicronsPerPixel / degreesPerPixel(2);
            
            % we need to create the scene XY grid in retinal microns, not
            % scene microns, because the sensor is specified in retinal microns
            sceneSpatialSupportInMicrons = sceneGet(obj.scenePrivate,'spatial support','microns');
            degreesPerSample = sceneGet(obj.scenePrivate,'deg per samp');
            micronsPerSample = sceneGet(scene,'distPerSamp','microns');
            sceneSpatialSupportInDegrees(:,:,1) = sceneSpatialSupportInMicrons(:,:,1) / micronsPerSample(1) * degreesPerSample;
            sceneSpatialSupportInDegrees(:,:,2) = sceneSpatialSupportInMicrons(:,:,2) / micronsPerSample(2) * degreesPerSample;
            
            % spatial support in retinal microns
            sceneSpatialSupport(:,:,1) = sceneSpatialSupportInDegrees(:,:,1) * retinalMicronsPerDegreeX;
            sceneSpatialSupport(:,:,2) = sceneSpatialSupportInDegrees(:,:,2) * retinalMicronsPerDegreeY;
            
            obj.sceneXgrid = squeeze(sceneSpatialSupport(:,:,1)); 
            obj.sceneYgrid = squeeze(sceneSpatialSupport(:,:,2));
            obj.sceneXdata = squeeze(obj.sceneXgrid(1,:));  % x-positions from 1st row
            obj.sceneYdata = squeeze(obj.sceneYgrid(:,1));  % y-positions from 1st col
        end
        
        function setCurrentTimeBin(obj, currentTimeBin)
            set(obj.timeSlider, 'value', currentTimeBin);
            obj.updateOpticalImageDisplay(currentTimeBin);
            obj.updateSensorViewDisplay(currentTimeBin);
            obj.updateOuterSegmentResponseDisplays(currentTimeBin);
        end
        
        function generateVideo(obj, frameStepInMilliseconds, videoFilename)
            writerObj = VideoWriter(videoFilename, 'MPEG-4'); % H264 format
            writerObj.FrameRate = 60; 
            writerObj.Quality = 100;
            writerObj.open();
            timeStepInSeconds = obj.outerSegmentResponseTimeData(2)-obj.outerSegmentResponseTimeData(1);
            step = round(frameStepInMilliseconds/(timeStepInSeconds*1000));
            for currentTimeBin = 1:step:numel(obj.outerSegmentResponseTimeData)
                obj.setCurrentTimeBin(currentTimeBin);
                writerObj.writeVideo(getframe(obj.hFig));
            end
            % close video stream and save movie
            writerObj.close();
            fprintf('Video exported in %s/%s.\n', pwd,videoFilename);
        end
        
        function resizeWindow(obj, desiredSize)   
            oldPosition = get(obj.hFig,'Position');
            
            if (numel(desiredSize) == 2)
                set(obj.hFig,'Position',[oldPosition(1) oldPosition(2) desiredSize(1) desiredSize(2)]);
                return;
            end
            
            if (strcmp(obj.figOrientation, 'horizontalLayout'))
                width = desiredSize;
                height = width/obj.aspectRatio;
            else
                height = desiredSize;
                width = height*obj.aspectRatio;
            end
            
            set(obj.hFig,'Position',[oldPosition(1) oldPosition(2) width height]);
        end
    end
        
    methods (Access = private)  
        
        function init(obj, figureNo, figureTitle, figOrientation)
            obj.hFig = figure(figureNo);
            obj.figOrientation = figOrientation;
            clf(obj.hFig);
            if (strcmp(obj.figOrientation, 'horizontalLayout'))
                obj.aspectRatio = 1024/500;
                screenSize = get(0,'ScreenSize');
                screenSize(3) = screenSize(3)*0.65;
            elseif (strcmp(obj.figOrientation, 'verticalLayout'))
                obj.aspectRatio = 800/1000;
                screenSize = get(0,'ScreenSize');
                screenSize(4) = screenSize(4)*0.9;
            else
                error('Unknown figure orientation: ''%s''.', figOrientation);
            end
            if (strcmp(obj.figOrientation, 'horizontalLayout'))
                set(obj.hFig, 'Name', figureTitle, 'Menubar', 'none', 'Toolbar', 'none', 'Color', [0.1 0.1 0.1], ...
                    'Position',[10+rand(1,1)*300 10+rand(1,1)*100 screenSize(3) screenSize(3)/obj.aspectRatio]);
            else
                set(obj.hFig, 'Name', figureTitle, 'Menubar', 'none', 'Toolbar', 'none', 'Color', [0.1 0.1 0.1], ...
                    'Position',[10+rand(1,1)*300 10+rand(1,1)*100 screenSize(4)*obj.aspectRatio screenSize(4)]);
            end
            %set(obj.hFig,'SizeChangedFcn', {@resizeOSwindow, obj, aspectRatio});
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
            ylabel(c,'current (pAmps)');
            colormap(obj.axesStruct.outerSegmentXTresponseAxes, 'bone');
            
            set(obj.axesStruct.outerSegmentXTresponseAxes, ...
                 'XLim', [obj.outerSegmentResponseTimeData(1) obj.outerSegmentResponseTimeData(end)], ...
                 'YLim', [obj.outerSegmentResponseXYpositionData(1) obj.outerSegmentResponseXYpositionData(end)], ...
                 'CLim', obj.outerSegmentDisplayedCurrentRange, ...
                 'XColor', [1 1 1], 'YColor', [1 1 1]);
            ylabel(obj.axesStruct.outerSegmentXTresponseAxes,  'cone #', 'FontSize', 10);
            set(obj.axesStruct.outerSegmentXTresponseAxes, 'FontSize', 10);

            % The traces plot
            cla(obj.axesStruct.outerSegmentTracesAxes);
            hold(obj.axesStruct.outerSegmentTracesAxes, 'on');
            
            % find the L,M, and S cone with the strongest response
            [lConeIndex, mConeIndex, sConeIndex, ...
             bestLconeSensorRowPosition, bestLconeSensorColPosition, ...
             bestMconeSensorRowPosition, bestMconeSensorColPosition, ...
             bestSconeSensorRowPosition, bestSconeSensorColPosition] = obj.findStrongestResponsingCones();
            
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
            xlabel(obj.axesStruct.outerSegmentTracesAxes, 'time (seconds)', 'FontSize', 10);
            ylabel(obj.axesStruct.outerSegmentTracesAxes, 'current (pAmps)', 'FontSize', 10);
            set(obj.axesStruct.outerSegmentTracesAxes, 'FontSize', 10);
            
            % The XY plot
            obj.computeSpatiallyInterpolatedOuterSegment2DResponseMap(positionIndex);

            cla(obj.axesStruct.outerSegmentXYresponseAxes);
            obj.outerSegmentOverlayPlots.p3 = ...
                 imagesc('XData', obj.outerSegmentResponseHiResXpositionData, ...
                         'YData', obj.outerSegmentResponseHiResYpositionData, ..., 
                         'CData', obj.outerSegmentSpatiallyInterpolated2DResponseMap, ...
                         'parent', obj.axesStruct.outerSegmentXYresponseAxes);

            hold(obj.axesStruct.outerSegmentXYresponseAxes, 'on');
            plot(obj.axesStruct.outerSegmentXYresponseAxes, obj.outerSegmentResponseXpositionData(bestLconeSensorColPosition), obj.outerSegmentResponseYpositionData(bestLconeSensorRowPosition), 'ro');
            plot(obj.axesStruct.outerSegmentXYresponseAxes, obj.outerSegmentResponseXpositionData(bestMconeSensorColPosition), obj.outerSegmentResponseYpositionData(bestMconeSensorRowPosition), 'go');
            plot(obj.axesStruct.outerSegmentXYresponseAxes, obj.outerSegmentResponseXpositionData(bestSconeSensorColPosition), obj.outerSegmentResponseYpositionData(bestSconeSensorRowPosition), 'bo');
            hold(obj.axesStruct.outerSegmentXYresponseAxes, 'off');
            set(obj.axesStruct.outerSegmentXYresponseAxes, ...
                'XLim', [0 max(obj.outerSegmentResponseHiResXpositionData)], ...
                'YLim', [0 max(obj.outerSegmentResponseHiResYpositionData)], ...
                'CLim', obj.outerSegmentDisplayedCurrentRange, ...
                'XColor', [1 1 1], 'YColor', [1 1 1]);
            
            axis(obj.axesStruct.outerSegmentXYresponseAxes,'ij'); axis(obj.axesStruct.outerSegmentXYresponseAxes,'equal');
            axis(obj.axesStruct.outerSegmentXYresponseAxes, 'off');
            set(obj.axesStruct.outerSegmentXYresponseAxes, 'FontSize', 10, 'XTick', [], 'YTick', []);
            colormap(obj.axesStruct.outerSegmentXYresponseAxes, bone(1024));
            title(obj.axesStruct.outerSegmentXYresponseAxes, sprintf('t = %2.3f msec', 1000.0*obj.outerSegmentResponseTimeData(positionIndex)), 'Color', [0.9 0.7 0.1], 'FontSize', 12);
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
            title(obj.axesStruct.outerSegmentXYresponseAxes, sprintf('t = %2.3f msec', 1000.0*obj.outerSegmentResponseTimeData(kPos)), 'Color', [0.9 0.7 0.1], 'FontSize', 12);
        end
        
        function initSensorStimulusEncodingDisplay(obj)
            positionIndex = 10;
            currentSensorPosition = squeeze(obj.sensorPositionsInMicrons(positionIndex,:));
            obj.findScenePixelsUnderSensor(currentSensorPosition);
        end
        
        function initSensorViewDisplay(obj)
            positionIndex = 10;
            currentSensorPosition = squeeze(obj.sensorPositionsInMicrons(positionIndex,:));
            
            if strcmp(obj.zoomedInView, 'optical image')
                obj.findOpticalImagePixelsUnderSensor(currentSensorPosition);
            elseif strcmp(obj.zoomedInView, 'scene')
                obj.findScenePixelsUnderSensor(currentSensorPosition);
            end

            cla(obj.axesStruct.sensorViewAxes);
            obj.sensorViewOverlayPlots.p1 = image('XData', obj.sensorViewXdata, 'YData', obj.sensorViewYdata, 'CData', obj.sensorView, 'parent', obj.axesStruct.sensorViewAxes);
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
            set(obj.axesStruct.sensorViewAxes, 'FontSize', 10);
            tickPositions = -2000:20:2000;
            set(obj.axesStruct.sensorViewAxes, 'XTick', tickPositions, 'YTick', tickPositions);
            box(obj.axesStruct.sensorViewAxes, 'on');
            axis(obj.axesStruct.sensorViewAxes, 'off');
        end
        
        function updateSensorViewDisplay(obj, kPos)
            currentSensorPosition = squeeze(obj.sensorPositionsInMicrons(kPos,:));
           
            if strcmp(obj.zoomedInView, 'optical image')
                obj.findOpticalImagePixelsUnderSensor(currentSensorPosition);
            elseif strcmp(obj.zoomedInView, 'scene')
                obj.findScenePixelsUnderSensor(currentSensorPosition);
            end
            
            set(obj.sensorViewOverlayPlots.p1, 'XData', obj.sensorViewXdata, 'YData', obj.sensorViewYdata, 'CData', obj.sensorView);  
            xpos = currentSensorPosition(1) -obj.sensorSizeInMicrons(1)/2 +  obj.sensorXsamplingGrid;
            ypos = currentSensorPosition(2) -obj.sensorSizeInMicrons(2)/2 +  obj.sensorYsamplingGrid;
            set(obj.sensorViewOverlayPlots.p2, 'XData', xpos, 'YData', ypos); 
            set(obj.axesStruct.sensorViewAxes, ...
                 'XLim', round(currentSensorPosition(1) + obj.sensorSizeInMicrons(1)*0.55*[-1 1]), ...
                 'YLim', round(currentSensorPosition(2) + obj.sensorSizeInMicrons(2)*0.55*[-1 1]));
        end
        
        function findScenePixelsUnderSensor(obj, currentSensorPosition)
            % find scene pixels falling within the sensor outline
            pixelIndices = find(...
                (obj.sceneXgrid >= currentSensorPosition(1) - obj.sensorSizeInMicrons(1)*0.6) & ...
                (obj.sceneXgrid <= currentSensorPosition(1) + obj.sensorSizeInMicrons(1)*0.6) & ...
                (obj.sceneYgrid >= currentSensorPosition(2) - obj.sensorSizeInMicrons(2)*0.6) & ...
                (obj.sceneYgrid <= currentSensorPosition(2) + obj.sensorSizeInMicrons(2)*0.6) );
            [rows, cols] = ind2sub(size(obj.sceneXgrid), pixelIndices);
            
            rowRange = min(rows):1:max(rows);
            colRange = min(cols):1:max(cols);
            obj.sensorView = obj.sceneRGBrenderingFullRes(rowRange,colRange,:);

            xGridSubset = obj.sceneXgrid(rowRange, colRange);
            yGridSubset = obj.sceneYgrid(rowRange, colRange);
            obj.sensorViewXdata = squeeze(xGridSubset(1,:));
            obj.sensorViewYdata = squeeze(yGridSubset(:,1));
        end
        
        function findOpticalImagePixelsUnderSensor(obj,currentSensorPosition)
            % find optical image pixels falling within the sensor outline
            pixelIndices = find(...
                (obj.opticalImageXgrid >= currentSensorPosition(1) - obj.sensorSizeInMicrons(1)*0.6) & ...
                (obj.opticalImageXgrid <= currentSensorPosition(1) + obj.sensorSizeInMicrons(1)*0.6) & ...
                (obj.opticalImageYgrid >= currentSensorPosition(2) - obj.sensorSizeInMicrons(2)*0.6) & ...
                (obj.opticalImageYgrid <= currentSensorPosition(2) + obj.sensorSizeInMicrons(2)*0.6) );
            [rows, cols] = ind2sub(size(obj.opticalImageXgrid), pixelIndices);
            
            rowRange = min(rows):1:max(rows);
            colRange = min(cols):1:max(cols);
            obj.sensorView = obj.opticalImageRGBrenderingFullRes(rowRange,colRange,:);
            xGridSubset = obj.opticalImageXgrid(rowRange, colRange);
            yGridSubset = obj.opticalImageYgrid(rowRange, colRange);
            obj.sensorViewXdata = squeeze(xGridSubset(1,:));
            obj.sensorViewYdata = squeeze(yGridSubset(:,1));
        end
        
        function [lConeIndex, mConeIndex, sConeIndex, ...
             bestLconeSensorRowPosition, bestLconeSensorColPosition, ...
             bestMconeSensorRowPosition, bestMconeSensorColPosition, ...
             bestSconeSensorRowPosition, bestSconeSensorColPosition] = findStrongestResponsingCones(obj)
         
            % get the cone types
            coneTypes = sensorGet(obj.sensorPrivate, 'cone type');
            % find the L,M, and Scone indices
            lConeIndices = find(coneTypes == 2);
            mConeIndices = find(coneTypes == 3);
            sConeIndices = find(coneTypes == 4);
            
            % find the strongest responding Lcone
            A = obj.outerSegmentXTCurrent(lConeIndices,:);
            [~,index] = max(A(:)); [row, ~] = ind2sub(size(A),index);
            lConeIndex = lConeIndices(row);
            
            % find the strongest responding Mcone
            A = obj.outerSegmentXTCurrent(mConeIndices,:);
            [~,index] = max(A(:)); [row, ~] = ind2sub(size(A),index);
            mConeIndex = mConeIndices(row);
            
            % find the strongest responding Scone
            A = obj.outerSegmentXTCurrent(sConeIndices,:);
            [~,index] = max(A(:)); [row, ~] = ind2sub(size(A),index);
            sConeIndex = sConeIndices(row);
            
            % Find the sensor (row,col) position for the strongest responding L,M and S-cones
            coneTypes = sensorGet(obj.sensorPrivate, 'cone type');
            [bestLconeSensorRowPosition, bestLconeSensorColPosition] = ind2sub(size(coneTypes), lConeIndex);
            [bestMconeSensorRowPosition, bestMconeSensorColPosition] = ind2sub(size(coneTypes), mConeIndex);
            [bestSconeSensorRowPosition, bestSconeSensorColPosition] = ind2sub(size(coneTypes), sConeIndex);
        end
        
        function computeSpatiallyInterpolatedOuterSegment2DResponseMap(obj, kPos)
            
            delta = obj.outerSegmentResponseHiResXpositionData(2)-obj.outerSegmentResponseHiResXpositionData(1); stepSize = round(1.0/delta);
            % get instant XY slice
            instantActivation = squeeze(obj.outerSegmentXYTCurrent(:,:,kPos));
            % normalize to [0 .. 1]
            minActivation = min(instantActivation(:))+0.01; maxActivation = max(instantActivation(:));
            instantActivation = (instantActivation-minActivation)/(maxActivation-minActivation);
            % fill zero padded array
            obj.outerSegmentSpatiallyInterpolated2DResponseMap = zeros(numel(obj.outerSegmentResponseHiResYpositionData), numel(obj.outerSegmentResponseHiResXpositionData));
            obj.outerSegmentSpatiallyInterpolated2DResponseMap(1 + (1:size(obj.outerSegmentXYTCurrent,1))*stepSize, 1+(1:size(obj.outerSegmentXYTCurrent,2))*stepSize) = instantActivation;
            % interpolate via convolution
            obj.outerSegmentSpatiallyInterpolated2DResponseMap = conv2(obj.outerSegmentSpatiallyInterpolated2DResponseMap, obj.outerSegmentXYResponseInterpolatingKernel, 'same');
            indices = find(abs(obj.outerSegmentSpatiallyInterpolated2DResponseMap - obj.outerSegmentSpatiallyInterpolated2DResponseMap(1,1))< 100*eps);
            % scale back to original scale
            obj.outerSegmentSpatiallyInterpolated2DResponseMap = obj.outerSegmentSpatiallyInterpolated2DResponseMap * (maxActivation-minActivation) + minActivation;
            obj.outerSegmentSpatiallyInterpolated2DResponseMap(indices) = obj.outerSegmentDisplayedCurrentRange(1)-1;
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
            sensorPositionHistory = obj.sensorPositionsInMicrons(max([1 kPos-1000]):kPos,:);
            set(obj.opticalImageOverlayPlots.p1, 'XData', currentSensorPosition(1) + obj.sensorOutlineInMicrons(:,1), 'YData', currentSensorPosition(2) + obj.sensorOutlineInMicrons(:,2));
            set(obj.opticalImageOverlayPlots.p2, 'XData', currentSensorPosition(1) + obj.sensorOutlineInMicrons(:,1), 'YData', currentSensorPosition(2) + obj.sensorOutlineInMicrons(:,2));
            set(obj.opticalImageOverlayPlots.p3, 'XData', sensorPositionHistory(:,1), 'YData', sensorPositionHistory(:,2));
        end
        
        function updateDisplaysWithNewDisplayedCurrentRange(obj)
            set(obj.axesStruct.outerSegmentXTresponseAxes, 'CLim', obj.outerSegmentDisplayedCurrentRange);
            set(obj.axesStruct.outerSegmentXYresponseAxes, 'CLim', obj.outerSegmentDisplayedCurrentRange);
            set(obj.axesStruct.outerSegmentTracesAxes, 'YLim', obj.outerSegmentDisplayedCurrentRange);
            set(obj.outerSegmentOverlayPlots.p7, 'YData', obj.outerSegmentDisplayedCurrentRange);
        end

        function generateAxesAndControls(obj, figOrientation) 
            if (strcmp(figOrientation, 'horizontalLayout'))
                w = 1024;
                h = 500;
            else
                w = 800;
                h = 1000;
            end
            imageWidthToHeightRatio = size(obj.opticalImageRGBrendering,2) / size(obj.opticalImageRGBrendering,1);
            
            leftMargin = 5/w;
            
            
            if (strcmp(figOrientation, 'horizontalLayout'))
                opticalImageWidth  = (w*0.6-10)/w;
                opticalImageHeight = (w*0.6-10)/imageWidthToHeightRatio/h;
                
                sensorViewWidth = 150/w; sensorViewHeight = 150/h; 
                spatiotemporalViewWidth = 380/w; spatiotemporalViewHeight = 150/h;
                spatialViewWidth  = 150/w; spatialViewHeight = 150/h; 
                
                xLeft = 7*leftMargin + opticalImageWidth;
                xLeft2 = xLeft;
                xLeft3 = xLeft2 + 130/w;
                bottomMargin1 = (h-10)/h - opticalImageHeight + 5/h;
                bottomMargin2 = (h-0)/h;
                bottomMargin3 = (h-10)/h - opticalImageHeight/2 -1.5*spatiotemporalViewHeight+120/h;
                bottomMargin4 = (h-10)/h - opticalImageHeight/2 -1.5*spatiotemporalViewHeight+15/h;
            else
                opticalImageWidth  = (w-10)/w;
                opticalImageHeight = (w-10)/imageWidthToHeightRatio/h;
                
                sensorViewWidth = 200/w; sensorViewHeight = 200/h; 
                spatiotemporalViewWidth = 500/w; spatiotemporalViewHeight = 200/h;
                spatialViewWidth  = 200/w; spatialViewHeight = 200/h; 
                
                xLeft = 9*leftMargin;
                xLeft2 = xLeft + 30/w + spatiotemporalViewWidth;
                xLeft3 = xLeft2;
                bottomMargin1 = (h-10)/h - opticalImageHeight + 5/h;
                bottomMargin2 = bottomMargin1;
                bottomMargin3 = bottomMargin1;
                bottomMargin4 = bottomMargin3-1.5*spatiotemporalViewHeight-25/h;
            end
            
            % generate plot axes
            obj.axesStruct.opticalImageAxes = axes('parent',obj.hFig,'unit','normalized','position',[leftMargin bottomMargin1 opticalImageWidth opticalImageHeight], 'Color', [0 0 0]);
            obj.axesStruct.sensorViewAxes   = axes('parent',obj.hFig,'unit','normalized','position',[leftMargin+20/w bottomMargin1+20/h sensorViewWidth sensorViewHeight], 'Color', [0 0 0]);

            
            % generate response axes
            obj.axesStruct.outerSegmentXTresponseAxes = axes('parent',obj.hFig, 'unit','normalized','position',[xLeft bottomMargin2-spatiotemporalViewHeight-10/h     spatiotemporalViewWidth spatiotemporalViewHeight], 'Color', [0 0 0]);
            obj.axesStruct.outerSegmentTracesAxes = axes('parent',    obj.hFig, 'unit','normalized','position',[xLeft bottomMargin2-1.5*spatiotemporalViewHeight-25/h spatiotemporalViewWidth spatiotemporalViewHeight/2], 'Color', [0 0 0]);
            
            % generate 2D instantaneous response axes
            positionVector = [xLeft2 bottomMargin4 spatialViewWidth spatialViewHeight];
            obj.axesStruct.outerSegmentXYresponseAxes = axes('parent',obj.hFig,'unit','normalized','position', positionVector, 'Color', [0 0 0]);
            
            % generate the optical image/scene textedit
            sensorViewColor = [0.7 0.7 0.6];
            %positionVector = [5*leftMargin+67/w+spatiotemporalViewWidth bottomMargin-0.5*spatialViewHeight+60/h 0.20 0.020];
            positionVector = [xLeft3+50/w bottomMargin3-0.5*spatialViewHeight+59/h 80/w 0.020];
            obj.sensorViewText = uicontrol(...
                'Parent', obj.hFig,...
                'String', sprintf('sensorView: opt.image'), ...
                'Style', 'edit',...
                'BackgroundColor', [0.1 0.1 0.1], ...
                'ForegroundColor',  sensorViewColor, ...
                'HorizontalAlignment', 'right', ...
                'FontSize', 10,...
                'Enable', 'inactive', ...
                'Units', 'normalized',...
                'Position', positionVector);
            
            % generate the optical image/scene slider
            positionVector = [xLeft3+135/w bottomMargin3-0.5*spatialViewHeight+57/h w/10000 0.02];
            obj.sensorViewSlider = uicontrol(...
                'Parent', obj.hFig,...
                'Style', 'slider',...
                'BackgroundColor',  sensorViewColor, ...
                'Min', 0, 'Max', 1, 'Value', 0, 'SliderStep', [1.0, 1.0], ...
                'Units', 'normalized',...
                'Position', positionVector);   
            % set the slider's callback function
            addlistener(obj.sensorViewSlider,'ContinuousValueChange', ...
                                      @(hFigure,eventdata) sensorViewSliderCallback(obj.sensorViewSlider,eventdata, obj));
            
            % generate the min and max displayed response editboxes
            displayedResponseRangeColor = [0.4 0.45 0.5]; 
            positionVector = [xLeft3+50/w bottomMargin3-0.5*spatialViewHeight+8/1000 80/w 0.020];
            uicontrol(...
                'Parent', obj.hFig,...
                'String', sprintf('max. pAmps'), ...
                'Style', 'edit',...
                'BackgroundColor', [0.1 0.1 0.1], ...
                'ForegroundColor', displayedResponseRangeColor, ...
                'HorizontalAlignment', 'right', ...
                'FontSize', 10,...
                'Enable', 'inactive', ...
                'Units', 'normalized',...
                'Position', positionVector);
            
            positionVector = [xLeft3+50/w bottomMargin3-0.5*spatialViewHeight+28/1000 80/w 0.020];
            uicontrol(...
                'Parent', obj.hFig,...
                'String', sprintf('min. pAmps'), ...
                'Style', 'edit',...
                'BackgroundColor', [0.1 0.1 0.1], ...
                'ForegroundColor',  displayedResponseRangeColor, ...
                'HorizontalAlignment', 'right', ...
                'FontSize', 10,...
                'Enable', 'inactive', ...
                'Units', 'normalized',...
                'Position', positionVector);
            
            positionVector = [xLeft3+135/w bottomMargin3-0.5*spatialViewHeight+25/1000 w/10000 0.02];
            obj.minDisplayedReponseSlider = uicontrol(...
                'Parent', obj.hFig,...
                'Style', 'slider',...
                'BackgroundColor',  displayedResponseRangeColor, ...
                'Min', -300, 'Max', 0, 'Value', -100,...
                'Units', 'normalized',...
                'Position', positionVector);   
            % set the slider's callback function
            addlistener(obj.minDisplayedReponseSlider,'ContinuousValueChange', ...
                                      @(hFigure,eventdata) minDisplayedResponseSliderCallback(obj.minDisplayedReponseSlider,eventdata, obj));
            
                                  
            positionVector = [xLeft3+135/w bottomMargin3-0.5*spatialViewHeight+5/1000 w/10000 0.02];
            obj.maxDisplayedReponseSlider = uicontrol(...
                'Parent', obj.hFig,...
                'Style', 'slider',...
                'BackgroundColor',  displayedResponseRangeColor, ...
                'Min', -80, 'Max', 100, 'Value', 0,...
                'Units', 'normalized',...
                'Position', positionVector);   
            % set the slider's callback function
            addlistener(obj.maxDisplayedReponseSlider,'ContinuousValueChange', ...
                                      @(hFigure,eventdata) maxDisplayedResponseSliderCallback(obj.maxDisplayedReponseSlider,eventdata, obj));
            
            % generate the time slider
            timeSliderLeftMargin = leftMargin;
            timeSliderBottom = (5)/h;
            
            obj.timeSlider = uicontrol(...
                'Parent', obj.hFig,...
                'Style', 'slider',...
                'BackgroundColor', [0.6 0.5 0.1], ...
                'Min', 1, 'Max', size(obj.sensorPositionsInMicrons,1), 'Value', 1,...
                'Units', 'normalized',...
                'Position', [timeSliderLeftMargin, timeSliderBottom 0.99 0.012]);    
           
            % set the slider's callback function
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
    Set(obj.hFig,'Position',[posVector(1:2) width height]);
end

% Callback for the sensor view slider
function sensorViewSliderCallback(hObject,eventdata, obj)
    if (get(hObject,'Value') < 0.5)
        set(obj.sensorViewText, 'String', 'sensorView: opt.image');
        obj.zoomedInView = 'optical image';
    else
        set(obj.sensorViewText, 'String', 'sensorView: scene');
        obj.zoomedInView = 'scene';
    end
    currentTimeBin = round(get(obj.timeSlider,'Value'));
    obj.updateSensorViewDisplay(currentTimeBin);
end

% Callback for time slider
function timeSliderCallback(hObject,eventdata, obj)
    currentTimeBin = round(get(hObject,'Value'));
    obj.updateOpticalImageDisplay(currentTimeBin);
    obj.updateSensorViewDisplay(currentTimeBin);
    obj.updateOuterSegmentResponseDisplays(currentTimeBin);
end


% Callback for minDisplayedResponseSlider
function minDisplayedResponseSliderCallback(hObject,eventdata, obj)
    newVal = round(get(hObject,'Value'));
    if (newVal <  obj.outerSegmentDisplayedCurrentRange(2))
        obj.outerSegmentDisplayedCurrentRange(1) = newVal;
        obj.updateDisplaysWithNewDisplayedCurrentRange();
    else
        fprintf(2, 'min displayed current cannot be larger than max displayed current...\n');
    end
end

% Callback for maxDisplayedResponseSlider
function maxDisplayedResponseSliderCallback(hObject,eventdata, obj)
    newVal = round(get(hObject,'Value'));
    if (newVal >  obj.outerSegmentDisplayedCurrentRange(1))
        obj.outerSegmentDisplayedCurrentRange(2) = newVal;
        obj.updateDisplaysWithNewDisplayedCurrentRange();
    else
        fprintf(2, 'max displayed current cannot be smaller than min displayed current...\n');
    end
end


