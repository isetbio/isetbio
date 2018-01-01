function [uData, hFig] = visualize(obj,varargin)
% Visualize an OI sequence
%
% Syntax
%   oiSequence.visualize;
%
% Description
%
% Inputs (required)
%
% Optional key/value pairs
%  format    - {'movie', 'weights', 'montage'} (default 'movie')
%  save      - Save a video file (boolean, false)
%  FrameRate - Frames per second (default 20)
%  vname  - Video file name for saving
%  showIlluminanceMap  -
%  eyeMovementsData    -  Show eye movement data (boolean, false)
%
% Return
%   uData - Displayed data
%   hFig  - Figure handle
%
% Parameter/value
%
% NP/BW ISETBIO Team, 2016

%% Interpret parameter values
p = inputParser;

p.addRequired('obj');

% For video case ...
p.addParameter('format','movieilluminance',@ischar);
p.addParameter('save',false,@islogical);
p.addParameter('vname','videoName',@ischar);
p.addParameter('FrameRate',20,@isnumeric);
p.addParameter('step',1,@isnumeric);
p.addParameter('showIlluminanceMap', false, @islogical);
p.addParameter('eyeMovementsData', struct('show', false), @(x)(isstruct(x)&&(isfield(x,'show'))));

varargin = ieParamFormat(varargin);
p.parse(obj,varargin{:});

format     = p.Results.format;
save       = p.Results.save;
vname      = p.Results.vname;
FrameRate  = p.Results.FrameRate;

%%  Show the oiSequence in one of the possible formats
uData = [];
vObj = [];    % Video object

switch format
    case 'weights'
        % Graph the weights'
        hFig = vcNewGraphWin;
        plot(obj.timeAxis, obj.modulationFunction);
        xlabel('Time (ms)'); ylabel('Modulation');
        title(sprintf('Composition: %s',obj.composition));
    case 'movieilluminance'
        % Show the oi as an illuminance movie
        wgts     = obj.modulationFunction;
        nFrames  = length(wgts);
        illFixed = oiGet(obj.oiFixed,'illuminance');
        illMod   = oiGet(obj.oiModulated,'illuminance');
        name     = oiGet(obj.oiModulated,'name');
        
        if save
            vObj = VideoWriter(vname);
            vObj.FrameRate = FrameRate;
            open(vObj);
        end
        
        % This code is general, and it could become an obj.get.movie;
        % Or obj.get.illuminanceMovie
        % The algorithm for mixing these is problematic because we
        % calculate the max between the two scenes.  This normalization can
        % lead to unwanted problems (as it did for vernier coding).  I need
        % to have the data come here in real physical units and deal with
        % it appropriately.
        mx1 = max(illFixed(:)); mx2 = max(illMod(:));
        mx = max(mx1,mx2);
        d = zeros([size(illFixed),length(obj.timeAxis)]);
        illFixed = 256*illFixed/mx; illMod = 256*illMod/mx;
        
        switch obj.composition
            case 'blend'
                for ii=1:nFrames
                    d(:,:,ii) = illFixed*(1-wgts(ii)) + illMod*wgts(ii);
                    % To make a video, we should do this type of thing
                end
            case 'add'
                for ii=1:nFrames
                    d(:,:,ii) = illFixed + illMod*wgts(ii);
                end     
            otherwise
                error('Unknown composition method: %s\n',obj.composition);
        end
        
        %  Show the movie data
        hFig = vcNewGraphWin; 
        colormap(gray(max(d(:)))); axis image; axis off;
        for ii=1:nFrames
            image(d(:,:,ii)); axis image; title(name); drawnow;
            if save,  F = getframe; writeVideo(vObj,F); end
        end
        
        % Write the video object if save is true
        if save
            writeVideo(vObj,F);
            close(vObj);
        end

        uData.movie = d;
    case 'moviergb'
        % Show the oi as an illuminance movie
        wgts     = obj.modulationFunction;
        nFrames  = length(wgts);
        rgbFixed = oiGet(obj.oiFixed,'rgb');
        rgbMod   = oiGet(obj.oiModulated,'rgb');
        name     = oiGet(obj.oiModulated,'name');
        
        if save
            vObj = VideoWriter(vname);
            vObj.FrameRate = FrameRate;
            open(vObj);
        end
        
        % This code is general, and it could become an obj.get.movie;
        % Or obj.get.illuminanceMovie
        % The algorithm for mixing these is problematic because we
        % calculate the max between the two scenes.  This normalization can
        % lead to unwanted problems (as it did for vernier coding).  I need
        % to have the data come here in real physical units and deal with
        % it appropriately.
        mx1 = max(rgbFixed(:)); mx2 = max(rgbMod(:));
        mx = max(mx1,mx2);
        d = zeros([size(rgbFixed),3,length(obj.timeAxis)]);
        rgbFixed = 256*rgbFixed/mx; illMod = 256*rgbMod/mx;
        
        switch obj.composition
            case 'blend'
                for ii=1:nFrames
                    d(:,:,:,ii) = rgbFixed*(1-wgts(ii)) + rgbMod*wgts(ii);
                    % To make a video, we should do this type of thing
                end
            case 'add'
                for ii=1:nFrames
                    d(:,:,:,ii) = rgbFixed + rgbMod*wgts(ii);
                end     
            otherwise
                error('Unknown composition method: %s\n',obj.composition);
        end
        
        %  Show the movie data
        hFig = vcNewGraphWin; 
        axis image; axis off;
        % Probably imagescRGB is not so good without a fixed scale!
        for ii=1:nFrames
            imagescRGB(d(:,:,:,ii)); axis image; title(name); drawnow;
            if save,  F = getframe; writeVideo(vObj,F); end
        end
        
        % Write the video object if save is true
        if save
            writeVideo(vObj,F);
            close(vObj);
        end

        uData.movie = d;
    case 'montage'
        % Window with snapshots
        colsNum = round(1.3*sqrt(obj.length));
        rowsNum = round(obj.length/colsNum);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum+1, ...
            'heightMargin',   0.05, ...
            'widthMargin',    0.02, ...
            'leftMargin',     0.04, ...
            'rightMargin',    0.00, ...
            'bottomMargin',   0.03, ...
            'topMargin',      0.03);
        
        if (p.Results.showIlluminanceMap)
            minIllum = Inf;
            maxIllum = -Inf;
            for oiIndex = 1:obj.length
                currentOI = obj.frameAtIndex(oiIndex);
                [illuminanceMap, ~] = oiCalculateIlluminance(currentOI);
                minIllum = min([minIllum min(illuminanceMap(:))]);
                maxIllum = max([maxIllum max(illuminanceMap(:))]);
            end
            if (minIllum == maxIllum)
                illumRange = [minIllum*0.99 maxIllum*1.01];
            else
                illumRange = [minIllum  maxIllum];
                meanIlluminance = mean(illumRange);
                illumMod = max(illumRange) / meanIlluminance - 1;
                illumRange = meanIlluminance + meanIlluminance*illumMod*[-1 1];
            end
        else
            XYZmax = 0;
            for oiIndex = 1:obj.length
                currentOI = obj.frameAtIndex(oiIndex);
                XYZ = oiGet(currentOI, 'xyz');
                if (max(XYZ(:)) > XYZmax)
                    XYZmax = max(XYZ(:));
                end
            end
            % Do not exceed XYZ values of 0.5 (for correct rendering)
            XYZmax = 2*XYZmax;
        end
        
        hFig = figure();
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 1700 730]); 
        for oiIndex = 1:obj.length
            if (oiIndex == 1)
                % Plot the modulation function
                subplot('Position', subplotPosVectors(1,1).v);
                bar(obj.timeAxis*1000, obj.modulationFunction, 0.9, 'LineWidth', 1.5, 'FaceColor', [1 0.5 0.5], 'EdgeColor', [1 0 0]);
                if (numel(obj.timeAxis)>1)
                    timeRange = [obj.timeAxis(1) obj.timeAxis(end)];
                else
                    timeRange = obj.timeAxis(1)+[-0.1 0.1];
                end
                set(gca, 'XLim', timeRange*1000, 'FontSize', 12);
                title(sprintf('composition: ''%s''', obj.composition));
                ylabel('modulation');
            end
            
            % Ask theOIsequence to return the oiIndex-th frame
            currentOI = obj.frameAtIndex(oiIndex);
            currentOIonsetTimeMillisecs = 1000*obj.timeAxis(oiIndex);
            dataXYZ = oiGet(currentOI, 'xyz');
            illuminanceMap = squeeze(dataXYZ(:,:,2));
            meanIlluminance = mean(illuminanceMap(:));
            %[illuminanceMap, meanIlluminance] = oiCalculateIlluminance(currentOI);
            support = oiGet(currentOI, 'spatial support', 'microns');
            xaxis = support(1,:,1);
            yaxis = support(:,1,2);
            row = 1+floor((oiIndex)/(colsNum+1));
            col = 1+mod((oiIndex),(colsNum+1));
            if (col > colsNum) || (row > rowsNum)
                continue;
            end
            subplot('Position', subplotPosVectors(row,col).v);
            if (p.Results.showIlluminanceMap)
                illuminanceMap = (illuminanceMap-illumRange(1))/(illumRange(2)-illumRange(1));
                imagesc(xaxis, yaxis, illuminanceMap);
                set(gca, 'CLim', [0 1]);
            else
                rgbImage = xyz2srgb(oiGet(currentOI, 'xyz')/XYZmax);
                imagesc(xaxis, yaxis, rgbImage, [0 1]);
            end

            axis 'image'
            if (col == 1) && (row == rowsNum)
                xticks = [xaxis(1) 0 xaxis(end)];
                yticks = [yaxis(1) 0 yaxis(end)];
                set(gca, 'XTick', xticks, 'YTick', yticks, 'XTickLabel', sprintf('%2.0f\n', xticks), 'YTickLabel', sprintf('%2.0f\n', yticks));
                ylabel('microns');
            else
                set(gca, 'XTick', [], 'YTick', [])
                xlabel(sprintf('frame %d (%2.1fms)', oiIndex, currentOIonsetTimeMillisecs));
            end
            
            if (p.Results.eyeMovementsData.show)
                hold on
                if (oiIndex < obj.length )
                    nextOIonsetTimeMillisecs = 1000*obj.timeAxis(oiIndex+1);
                else
                    nextOIonsetTimeMillisecs = 1000*(obj.timeAxis(oiIndex) +(obj.timeAxis(oiIndex)-obj.timeAxis(oiIndex-1)));
                end
            
                % plot eye movements during previous OIs in black
                idx = find(p.Results.eyeMovementsData.timeAxisMillisecs < currentOIonsetTimeMillisecs);
                plot(p.Results.eyeMovementsData.posMicrons(idx,1), p.Results.eyeMovementsData.posMicrons(idx,2), 'k.-');
                 % plot eye movements during current OI in red
                idx = find(...
                    p.Results.eyeMovementsData.timeAxisMillisecs >= currentOIonsetTimeMillisecs & ...
                    p.Results.eyeMovementsData.timeAxisMillisecs < nextOIonsetTimeMillisecs ...
                );
                plot(p.Results.eyeMovementsData.posMicrons(idx,1), p.Results.eyeMovementsData.posMicrons(idx,2), 'r.-');
                hold off;
            end
            
            if (p.Results.showIlluminanceMap)
                colormap(jet(1024));
            end

            title(sprintf('mean illum: %2.4f td', meanIlluminance));
            set(gca, 'FontSize', 12);
        end

    otherwise
        error('Unknown format %s\n',format);
end

end

