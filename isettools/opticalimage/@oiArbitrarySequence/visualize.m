function [uData, vObj] = visualize(obj, plotType, varargin)
% Visualize an OI sequence
%
% Syntax:
%   [uData, vObj] = visualize(obj, plotType, [varargin]);
%
% Description
%    This is a plot method for the oiSequence class.
%
% Inputs:
%    obj       - Object. The oi sequence object
%    plotType  - String. The type of plot. Options are:
%        'movie illuminance': Gray scale (luminance) video of the stimulus
%        'movie rgb':         Video of the stimuli
%        'weights':           Plot showing time series of weights
%        'montage':           A large montage of the frames and first panel
%                             of weights
%
% Outputs:
%    uData     - Struct. A structure containing the displayed data in a
%                matrix named movie.
%    vObj      - Object. Video object for the movie case.
%
% Optional key/value pairs:
%    save      - Boolean. Save the movie. Default false.
%    FrameRate - Numeric. Frames per second. Default 20.
%    vname     - String. Video file name when saving. Default 'videoName'.
%    showIlluminanceMap
%              - Boolean. Show illuminance map. Default false.
%    eyeMovementsData
%              -  Boolean. Show eye movement data. Default false.
%
% See Also:
%   t_oisCreate
%

% History:
%    06/28/2020 NPC, ISETBio team, 2020

%% Interpret parameter values
p = inputParser;
plotType = ieParamFormat(plotType);
varargin = ieParamFormat(varargin);

p.addRequired('obj');
validTypes = {'movieilluminance','moviergb','montage'};
p.addRequired('plottype', @(x)(ismember(x,validTypes)));

% For video case ...
p.addParameter('vname', '', @ischar);
p.addParameter('framerate', 20, @isnumeric);

% Must ask NP more about this
p.addParameter('showilluminancemap', false, @islogical);
p.addParameter('eyemovementsdata', struct('show', false), ...
    @(x)(isstruct(x)&&(isfield(x, 'show'))));

% Whether to use vcGraphWin or matlab's figure for rendering
p.addParameter('backendrenderer', 'vcGraphWin', ...
    @(x)(ischar(x)&&(ismember(x, {'vcGraphWin', 'figure'}))));

p.parse(obj, plotType, varargin{:});

vname = p.Results.vname;
FrameRate = p.Results.framerate;

save = false;
if ~isempty(vname), save = true; end

%%  Show the oiSequence in one of the possible formats
uData = [];  % Returned data.
vObj = [];   % Video object

nFrames = length(obj.timeAxis);
switch ieParamFormat(plotType)

    case 'movieilluminance'
        % Show the oi as an illuminance movie
        
        name = oiGet(obj.frameAtIndex(1), 'name');
        for ii = 1:nFrames
            d(:, :, ii) = oiGet(obj.frameAtIndex(ii), 'illuminance');
        end
                
        %  Show the movie data. 20Hz Frame rate.
        h = vcNewGraphWin;
        colormap(gray(max(d(:))));
        axis image;
        axis off;
        for ii = 1:nFrames
            image(d(:, :, ii));
            axis image;
            title(name);
            drawnow;
            pause(0.05);
        end
        delete(h);

        uData.movie = d;

        % Write the video object if save is true
        if save
            %disp('Saving video ...')
            [~, vObj] = ieMovie(uData.movie, ...
                'vname', vname, ...
                'FrameRate', FrameRate, ...
                'show', false);
            %disp('Done')
        end

    case 'moviergb'
        % Show the oi as an RGB movie
        for ii = 1:nFrames
            d(:, :, :, ii) = xyz2rgb(oiGet(obj.frameAtIndex(ii), 'xyz'));
        end
        d = d / max(d(:));   
        name = oiGet(obj.frameAtIndex(1), 'name');

        %  Show the movie data. 20Hz Frame rate.
        h = vcNewGraphWin;
        axis image;
        axis off;
        for ii = 1:nFrames
            % imagesc(d(:, :, :, ii), [0 256]);
            image(d(:, :, :, ii));
            axis image;
            title(name);
            drawnow;
            pause(0.05);
        end
        delete(h);

        uData.movie = d;

        % Write the video object if save is true
        if save
            %disp('Saving video')
            [~, vObj] = ieMovie(uData.movie, ...
                'vname', vname, ...
                'FrameRate', FrameRate, ...
                'show', false);
            %disp('Done')
        end

    case 'montage'
        % Window with snapshots and possibly eye movements.
        colsNum = round(1.3 * sqrt(obj.length));
        rowsNum = ceil(obj.length / colsNum);
        subplotPosVectors = NicePlot.getSubPlotPosVectors(...
            'rowsNum', rowsNum, ...
            'colsNum', colsNum+1, ...
            'heightMargin', 0.05, ...
            'widthMargin', 0.02, ...
            'leftMargin', 0.04, ...
            'rightMargin', 0.00, ...
            'bottomMargin', 0.03, ...
            'topMargin', 0.03);

        if (p.Results.showilluminancemap)
            minIllum = Inf;
            maxIllum = -Inf;
            for oiIndex = 1:obj.length
                currentOI = obj.frameAtIndex(oiIndex);
                [illuminanceMap, ~] = oiCalculateIlluminance(currentOI);
                minIllum = min([minIllum min(illuminanceMap(:))]);
                maxIllum = max([maxIllum max(illuminanceMap(:))]);
            end
            if (minIllum == maxIllum)
                illumRange = [minIllum * 0.99 maxIllum * 1.01];
            else
                illumRange = [minIllum  maxIllum];
                meanIlluminance = mean(illumRange);
                illumMod = max(illumRange) / meanIlluminance - 1;
                illumRange = meanIlluminance + ...
                    meanIlluminance * illumMod * [-1 1];
            end
        else
            XYZmax = 0;
            for oiIndex = 1:obj.length
                currentOI = obj.frameAtIndex(oiIndex);
                XYZ = oiGet(currentOI, 'xyz');
                if (max(XYZ(:)) > XYZmax), XYZmax = max(XYZ(:)); end
            end
            % Do not exceed XYZ values of 0.5 (for correct rendering)
            XYZmax = 2 * XYZmax;
        end

        if strcmp(p.Results.backendrenderer, 'vcGraphWin')
            h = vcNewGraphWin;
        else
            h = figure(1000); clf;
            uData.figHandle = h;
        end
        set(h, 'Color', [1 1 1], 'Position', [10 10 1700 730]);
        for oiIndex = 1:obj.length

            % Ask theOIsequence to return the oiIndex-th frame
            currentOI = obj.frameAtIndex(oiIndex);
            currentOIonsetTimeMillisecs = 1000 * obj.timeAxis(oiIndex);
            dataXYZ = oiGet(currentOI, 'xyz');
            illuminanceMap = squeeze(dataXYZ(:, :, 2));
            meanIlluminance = mean(illuminanceMap(:));
            % [illuminanceMap, meanIlluminance] = ...
            %    oiCalculateIlluminance(currentOI);
            support = oiGet(currentOI, 'spatial support', 'microns');
            xaxis = support(1, :, 1);
            yaxis = support(:, 1, 2);
            if (p.Results.eyemovementsdata.show)
                spatialRange = [-1 1] * 1.05 * ...
                    max([max(abs(xaxis(:))), max(abs(yaxis(:))), ...
                    max(abs(p.Results.eyemovementsdata.posMicrons))]);
            else
                spatialRange = [-1 1] * 1.05 * ...
                    max([max(abs(xaxis(:))) max(abs(yaxis(:)))]);
            end
            row = floor(oiIndex / (colsNum + 1)) + 1;
            col = mod(oiIndex, colsNum + 1) + 1;
            ax = subplot('Position', subplotPosVectors(row, col).v);
            
            if (p.Results.showilluminancemap)
                illuminanceMap = (illuminanceMap - illumRange(1)) / ...
                    (illumRange(2) - illumRange(1));
                imagesc(ax, xaxis, yaxis, illuminanceMap);
                set(ax, 'CLim', [0 1]);
            else
                rgbImage = xyz2srgb(oiGet(currentOI, 'xyz') / XYZmax);
                imagesc(ax,xaxis, yaxis, rgbImage, [0 1]);
            end

            axis(ax, 'image')
            set(ax, 'XTick', [], 'YTick', [])
            xlabel(ax,sprintf('frame %d (%2.1fms)', oiIndex, currentOIonsetTimeMillisecs));
           
            
            if (p.Results.eyemovementsdata.show)
                hold(ax, 'on')
                if (oiIndex < obj.length )
                    nextOIonsetTimeMillisecs = ...
                        1000 * obj.timeAxis(oiIndex+1);
                else
                    nextOIonsetTimeMillisecs = 1000 * ...
                        (obj.timeAxis(oiIndex) + ...
                        (obj.timeAxis(oiIndex) - obj.timeAxis(oiIndex-1)));
                end

                % plot eye movements during previous OIs in black
                idx = find(p.Results.eyemovementsdata.timeAxisMillisecs ...
                    < currentOIonsetTimeMillisecs);
                plot(ax, p.Results.eyemovementsdata.posMicrons(idx, 1), ...
                    p.Results.eyemovementsdata.posMicrons(idx, 2), 'k.-');
                 % plot eye movements during current OI in red
                idx = find(...
                    p.Results.eyemovementsdata.timeAxisMillisecs >= ...
                    currentOIonsetTimeMillisecs & ...
                    p.Results.eyemovementsdata.timeAxisMillisecs < ...
                    nextOIonsetTimeMillisecs);
                plot(ax, p.Results.eyemovementsdata.posMicrons(idx, 1), ...
                    p.Results.eyemovementsdata.posMicrons(idx, 2), 'r.-');
                hold(ax, 'off');
            end

            if (p.Results.showilluminancemap), colormap(ax, bone(1024)); end

            set(ax, 'XLim', spatialRange, 'YLim', spatialRange);
            axis(ax, 'xy')

            title(ax, sprintf('mean illum: %2.4f td', meanIlluminance));
            set(ax, 'FontSize', 12);
            drawnow;
        end

    otherwise
        error('Unknown plot type %s\n', plotType);
end

end