function [uData, hFig] = visualize(obj,varargin)
% Visualize aspects of the OI sequence
%
% Parameter/value
%  format - {weights, movie, montage}
%  save   - Save a video file
%
% NP/BW ISETBIO Team, 2016

%% Interpret parameter values
p = inputParser;

p.addRequired('obj');

% For video case ...
p.addParameter('format','movie',@ischar);
p.addParameter('save',false,@islogical);
p.addParameter('vname','videoName',@ischar);
p.addParameter('FrameRate',20,@isnumeric);
p.addParameter('step',1,@isnumeric);

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
    case 'movie'
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

    case 'montage'
        % Window with snapshots
        colsNum = round(1.3*sqrt(obj.length));
        rowsNum = ceil(obj.length/colsNum);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum+1, ...
            'heightMargin',   0.05, ...
            'widthMargin',    0.02, ...
            'leftMargin',     0.04, ...
            'rightMargin',    0.00, ...
            'bottomMargin',   0.03, ...
            'topMargin',      0.03);
        
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
        
        hFig = figure();
        set(hFig, 'Color', [1 1 1], 'Position', [10 10 1700 730]);
        
        for oiIndex = 1:obj.length
            if (oiIndex == 1)
                % Plot the modulation function
                subplot('Position', subplotPosVectors(1,1).v);
                stairs(1:obj.length, obj.modulationFunction, 'r', 'LineWidth', 1.5);
                set(gca, 'XLim', [1 obj.length], 'FontSize', 12);
                title(sprintf('composition: ''%s''', obj.composition));
                xlabel('frame index');
                ylabel('modulation');
            end
            
            % Ask theOIsequence to return the oiIndex-th frame
            currentOI = obj.frameAtIndex(oiIndex);
            support = oiGet(currentOI, 'spatial support', 'microns');
            [~, meanIlluminance] = oiCalculateIlluminance(currentOI);
            xaxis = support(1,:,1);
            yaxis = support(:,1,2);
            row = 1+floor((oiIndex)/(colsNum+1));
            col = 1+mod((oiIndex),(colsNum+1));
            
            subplot('Position', subplotPosVectors(row,col).v);
            rgbImage = xyz2srgb(oiGet(currentOI, 'xyz')/XYZmax);
            imagesc(xaxis, yaxis, rgbImage, [0 1]);
            axis 'image'
            if (col == 1) && (row == rowsNum)
                xticks = [xaxis(1) 0 xaxis(end)];
                yticks = [yaxis(1) 0 yaxis(end)];
                set(gca, 'XTick', xticks, 'YTick', yticks, 'XTickLabel', sprintf('%2.0f\n', xticks), 'YTickLabel', sprintf('%2.0f\n', yticks));
                ylabel('microns');
            else
                set(gca, 'XTick', [], 'YTick', [])
                xlabel(sprintf('frame %d (%2.1fms)', oiIndex, 1000*obj.timeAxis(oiIndex)));
            end
            title(sprintf('mean illum: %2.1f', meanIlluminance));
            set(gca, 'FontSize', 12);
        end
    otherwise
        error('Unknown format %s\n',format);
end

end

