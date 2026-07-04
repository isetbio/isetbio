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
%    xx/xx/16  NP/BW  ISETBIO Team, 2016
%    01/06/18  dhb    Don't print text into command window.
%                     If we want text sometimes, add key/value 'verbose'
%                     pair and set default to false.
%    03/27/18  jnm    Formatting (Also 07/01/19)

%% Interpret parameter values
p = inputParser;
plotType = ieParamFormat(plotType);
varargin = ieParamFormat(varargin);

p.addRequired('obj');
validTypes = {'movieilluminance','moviergb','weights','montage'};
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

switch ieParamFormat(plotType)
    case 'weights'
        % Graph the weights'
        vcNewGraphWin;
        plot(obj.timeAxis, obj.modulationFunction);
        xlabel('Time (ms)');
        ylabel('Weight');
        title(sprintf('Composition: %s', obj.composition));
        grid on;
        uData.time = obj.timeAxis;
        uData.wgts = obj.modulationFunction;
    case 'movieilluminance'
        % Show the oi as an illuminance movie
        wgts = obj.modulationFunction;
        nFrames = length(wgts);
        illFixed = oiGet(localEnsurePhotons(obj.oiFixed,    obj.photonsFixed),    'illuminance');
        illMod   = oiGet(localEnsurePhotons(obj.oiModulated, obj.photonsModulated), 'illuminance');
        name = oiGet(obj.oiModulated, 'name');

        mx1 = max(illFixed(:));
        mx2 = max(illMod(:));
        mx  = max(mx1, mx2);
        illFixed = 256 * illFixed / mx;
        illMod   = 256 * illMod   / mx;

        % Determine colormap scale without materializing the full movie.
        switch obj.composition
            case 'blend',  cmapMax = 256;
            case 'add',    cmapMax = ceil(256 * (1 + max(abs(wgts))));
            otherwise,     cmapMax = 256;
        end

        if save, d = zeros([size(illFixed) nFrames]); end

        h = ieFigure;
        colormap(gray(cmapMax));
        axis image; axis off;
        for ii = 1:nFrames
            switch obj.composition
                case 'blend'
                    frame = illFixed * (1 - wgts(ii)) + illMod * wgts(ii);
                case 'add'
                    frame = illFixed + illMod * wgts(ii);
                otherwise
                    error('Unknown composition method: %s\n', obj.composition);
            end
            if save, d(:, :, ii) = frame; end
            image(frame); axis image; title(name); drawnow; pause(0.05);
        end
        delete(h);

        if save
            uData.movie = d;
            [~, vObj] = ieMovie(d, 'vname', vname, 'FrameRate', FrameRate, ...
                'show', false);
        end

    case 'moviergb'
        % Show the oi as an RGB movie
        wgts = obj.modulationFunction;
        nFrames = length(wgts);

        % I am not sure why this does not work as well with
        % oiGet(oi, 'rgb');  There appears to be some scaling in that case
        % that shifts the means.
        xyzMod   = oiGet(localEnsurePhotons(obj.oiModulated, obj.photonsModulated), 'xyz');
        xyzFixed = oiGet(localEnsurePhotons(obj.oiFixed,    obj.photonsFixed),    'xyz');
        rgbMod = xyz2rgb(xyzMod);
        rgbFixed = xyz2rgb(xyzFixed);
        name = oiGet(obj.oiModulated, 'name');

        % Scale the RGB data to [0, 1] with a common scale factor
        mx1 = max(rgbFixed(:));
        mx2 = max(rgbMod(:));
        mx  = max(mx1, mx2);
        rgbFixed = rgbFixed / mx;
        rgbMod   = rgbMod   / mx;

        if save, d = zeros([size(rgbFixed) nFrames]); end

        h = vcNewGraphWin;
        axis image; axis off;
        for ii = 1:nFrames
            switch obj.composition
                case 'blend'
                    % Means of rgbFixed and rgbMod should agree in this case.
                    frame = rgbFixed * (1 - wgts(ii)) + rgbMod * wgts(ii);
                case 'add'
                    frame = rgbFixed + rgbMod * wgts(ii);
                otherwise
                    error('Unknown composition method: %s\n', obj.composition);
            end
            if save, d(:, :, :, ii) = frame; end
            image(frame); axis image; title(name); drawnow; pause(0.05);
        end
        delete(h);

        if save
            uData.movie = d;
            [~, vObj] = ieMovie(d, 'vname', vname, 'FrameRate', FrameRate, ...
                'show', false);
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

        % Pass 1: find normalisation. photonsFixed/photonsModulated are
        % cached after the first frameAtIndex call, so pass 2 is fast.
        if (p.Results.showilluminancemap)
            minIllum = Inf;
            maxIllum = -Inf;
            for oiIndex = 1:obj.length
                currentOI = obj.frameAtIndex(oiIndex);
                [illuminanceMap, ~] = oiCalculateIlluminance(currentOI);
                minIllum = min(minIllum, min(illuminanceMap(:)));
                maxIllum = max(maxIllum, max(illuminanceMap(:)));
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
            h = figure();
            uData.figHandle = h;
        end
        set(h, 'Color', [1 1 1], 'Units', 'normalized', 'Position', [0.01 0.05 0.97 0.88]);
        figure(h);  % ensure h is the active figure before any subplot calls
        for oiIndex = 1:obj.length
            if (oiIndex == 1)
                % Plot the modulation function
                subplot('Position', subplotPosVectors(1, 1).v);
                bar(obj.timeAxis * 1000, obj.modulationFunction, 0.9, ...
                    'LineWidth', 1.5, 'FaceColor', [1 0.5 0.5], ...
                    'EdgeColor', [1 0 0]);
                if (numel(obj.timeAxis)>1)
                    timeRange = [obj.timeAxis(1) obj.timeAxis(end)];
                else
                    timeRange = obj.timeAxis(1) + [-0.1 0.1];
                end
                set(gca, 'XLim', timeRange * 1000, 'FontSize', 12);
                title(sprintf('composition: ''%s''', obj.composition));
                ylabel('modulation');
            end

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
            subplot('Position', subplotPosVectors(row, col).v);
            if (p.Results.showilluminancemap)
                illuminanceMap = (illuminanceMap - illumRange(1)) / ...
                    (illumRange(2) - illumRange(1));
                imagesc(xaxis, yaxis, illuminanceMap);
                set(gca, 'CLim', [0 1]);
            else
                rgbImage = xyz2srgb(oiGet(currentOI, 'xyz') / XYZmax);
                imagesc(xaxis, yaxis, rgbImage, [0 1]);
            end

            axis 'image'
            if (col == 1) && (row == rowsNum)
                xticks = [xaxis(1) 0 xaxis(end)];
                yticks = [yaxis(1) 0 yaxis(end)];
                set(gca, 'XTick', xticks, 'YTick', yticks, ...
                    'XTickLabel', sprintf('%2.0f\n', xticks), ...
                    'YTickLabel', sprintf('%2.0f\n', yticks));
                ylabel('microns');
            else
                set(gca, 'XTick', [], 'YTick', [])
                xlabel(sprintf('frame %d (%2.1fms)', oiIndex, ...
                    currentOIonsetTimeMillisecs));
            end

            if (p.Results.eyemovementsdata.show)
                hold on
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
                plot(p.Results.eyemovementsdata.posMicrons(idx, 1), ...
                    p.Results.eyemovementsdata.posMicrons(idx, 2), 'k.-');
                 % plot eye movements during current OI in red
                idx = find(...
                    p.Results.eyemovementsdata.timeAxisMillisecs >= ...
                    currentOIonsetTimeMillisecs & ...
                    p.Results.eyemovementsdata.timeAxisMillisecs < ...
                    nextOIonsetTimeMillisecs);
                plot(p.Results.eyemovementsdata.posMicrons(idx, 1), ...
                    p.Results.eyemovementsdata.posMicrons(idx, 2), 'r.-');
                hold off;
            end

            if (p.Results.showilluminancemap), colormap(bone(1024)); end

            set(gca, 'XLim', spatialRange, 'YLim', spatialRange);
            axis 'xy'

            title(sprintf('mean illum: %2.4f td', meanIlluminance));
            set(gca, 'FontSize', 12);
        end

    otherwise
        error('Unknown plot type %s\n', plotType);
end

end

% -------------------------------------------------------------------------

function oi = localEnsurePhotons(oi, photonsCache)
% Return oi with photons present.
% After the first frameAtIndex call, photons are moved from oiFixed/
% oiModulated into the cache (photonsFixed/photonsModulated) to eliminate
% the in-memory duplicate.  This helper restores them temporarily for any
% oiGet call that requires photons (illuminance, xyz, etc.).
if ~isempty(photonsCache)
    oi = oiSet(oi, 'photons', photonsCache);
end
end
