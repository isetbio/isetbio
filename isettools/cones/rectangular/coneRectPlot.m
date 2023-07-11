function [uData, hf] = coneRectPlot(obj, plotType, varargin)
% Plot function for @coneMosaicRect class
%
% Syntax:
%   [uData, hf] = plot(obj, plotType, app, varargin)
%
% Description:
%    There is a specialized plot method for the coneMosaic class that
%    calls this function.
%
%    When the plot type string begins with 'os ' or 'outersegment ' we pass
%    the arguments along to os.plot(). For example,
%
%           cMosaic.plot('os impulse response')
%
%    graphs the outer segment impulse response on its own time axis.
%
% Inputs:
%    obj      - Either the coneMosaicRect or the coneRectWindow_app,
%               which contains the coneMosaicRect in the cMosaic
%               slot
%    plotType - string, type of plot (N.B. The '*' means user must select a
%               point on an image to specify the location of the point or
%               line to be plotted.)
%       Cone mosaic              - Color image of the cone arrangement
%
%       Mean absorptions         - Image of the mean absorptions
%       Movie absorptions        - Movie of the absorptions on cone mosaic
%       Mean current             - Image of the mean current
%       Movie current            - Gray scale movie of current
%
%       hline absorptions*       - Graph of horizontal line of absoprtions
%       hline current*           -
%       hline absorptions lms*   - Three panel graph of LMS absorptions
%       hline current lms*       - Three panel graph of LMS current
%       vline absorptions*       - Vertical line
%       vline current*
%       vline absorptions lms*
%       vline current lms*
%       time series absorptions* - Cone absorptions Graph of a point
%       time series current*     - Cone photocurrent Graph of a point
%
%       Impulse response         - Cone current impulse response
%
%       Cone fundamentals        - Cone pigment without macular or lens
%       Cone spectral QE         - Cone pigment and macular
%       Eye spectral QE          - Cone pigment with macular and lens
%       Macular transmittance    - Graph
%       Macular absorptance      - Graph
%       Macular absorbance       - Graph
%       Eye movement path        - eye movement path
%
% Outputs:
%    uData    - The user data
%    hf       - The figure or axes handle
%
% Optional key/value pairs:
%    'hf' - figure handle or control structure, the meaning of value is
%           []: create plot in new figure
%           'none': don't plot
%         (isgraphics figure or axes handle): Use that figure or axes
%    'oi' - Optical image used to get lens transmittance for QE plots
%    'x' - x value (col) for vline plots
%    'y' - y value (row) for hline plots
%    roi - a point or rect for drawing a line
%
% Notes:
%    * TODO: Think about coneImageActivity function at end.
%

%  Examples:
%{
  cm = coneMosaicRect;
  coneRectPlot('help');
  coneRectPlot(cm,'cone mosaic');
%}
%{
  coneRectPlot(thisM,'hline absorptions','roi',[150 150]);
%}
%% If user wants help

%
validPlotsPrint = {'help', ...
    'Cone mosaic', ...
    'Mean absorptions', 'Movie absorptions', ...
    'Mean current', 'Movie current', ...
    'hline absorptions', 'hline current', ...
    'hline absorptions lms', 'hline current lms', ...
    'vline absorptions', 'vline current', ...
    'vline absorptions lms', 'vline current lms' ...
    'Impulse response', ...
    'time series absorptions', 'time series current', ...
    'Cone fundamentals', 'Cone spectral QE', 'Eye spectral QE', ...
    'Macular transmittance', 'Macular absorptance', ...
    'Macular absorbance', 'Eye movement path'};

% Squeeze spaces and force lower case on validPlots. Then check plotType
validPlots = cellfun(@(x)(ieParamFormat(x)), validPlotsPrint, ...
    'UniformOutput', false);

if ischar(obj) && strcmp(obj,'help')
    coneRectPlotHelp(validPlotsPrint);
    return;
end

%% Start parsing
varargin = ieParamFormat(varargin);
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('obj',...
    @(x)(isa(x,'coneMosaicRect') || ...
    isa(x,'coneRectWindow_App') || ...
    strcmp(obj,'help')));

p.addRequired('plotType',...
    @(x)(ismember(ieParamFormat(x),validPlots)) || ...
    @(x)(strcmpi(ieParamFormat(plotType(1:12)),'outersegment')));

p.addParameter('hf',[],@(x)(isa(x,'matlab.ui.Figure')));
p.addParameter('x', [], @isscalar);   % x axis value
p.addParameter('y', [], @isscalar);   % y axis value
p.addParameter('roi',[],@isvector);  % (x,y) or (x,y,width,height)
p.addParameter('oi',[],@isstruct);

p.parse(obj,plotType,varargin{:});

plotType = ieParamFormat(plotType);
hf  = p.Results.hf;
oi  = p.Results.oi;                    % Used in plotGraphs routine
roi = p.Results.roi;



%% Initialize the window for plotting

if isa(obj,'coneRectWindow_App')
    % User sent coneRectWindow_App.  Use the image for plotting,
    % unless the user sent in a handle to a figure, hf.
    app   = obj;   % Just renaming
    cm    = app.cMosaic;
    curAx = app.imgMain;
    if ~isempty(hf)
        curAx = get(hf,'CurrentAxes');
    end
else
    % User sent coneMosaicRect in
    app = []; cm  = obj;
    if isempty(hf)
        % If they did not specify a window, create one.
        hf = ieNewGraphWin;
    end
    curAx = get(hf,'CurrentAxes');
end


% Deal with the gamma from the window
if isa(app,'coneRectWindow_App'), gam = str2double(app.editGam.Value);
else, gam = 1;
end

% Check plot type to see if we send this off to the outer segment plot
% routine
if numel(plotType)> 11 && strcmpi(ieParamFormat(plotType(1:12)),'outersegment')
    cm.os.plot(plotType(13:end), 'cmosaic', obj, varargin{:});
end

%% Set color order so that LMS plots as RGB

% If the user has changed the default, we leave it alone.
% We should have a better way to check this.
co = get(curAx, 'ColorOrder');
if size(co,1) == 7    % This is the default number.  So we reorder.
    set(get(curAx, 'parent'), ...
        'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :))
end

%% Switch on passed plot type

% When we draw into the main axis in the cMosaic window, we store the user
% data in that window.  Otherwise, we build up some user data in a new
% variable that is stored in the temporary plotting window.
switch ieParamFormat(plotType)
    case 'conemosaic'
        % This might become
        %
        %   cmrPlot(obj,cone mosaic, axisData)
        %
        axisData = get(curAx,'UserData');  % Main axis in the window

        if isfield(axisData,'mosaicImage') && ~isempty(axisData.mosaicImage)
            imagesc(axisData.mosaicImage)
            axis off; axis image;
        else
            locs    = cm.coneLocs;
            pattern = cm.pattern(:);

            % The locations are converted to microns from meters, I think.
            [axisData.support, axisData.spread, axisData.delta, axisData.mosaicImage] = ...
                conePlot(locs * 1e6, pattern);
            imagesc(axisData.mosaicImage);
            axis off; axis image;
        end

        set(curAx,'UserData',axisData);  % Put the modified values back

    case 'meanabsorptions'
        % This might become
        %
        %   cmrPlot(obj,mean absorptions,axisData)
        %

        % Image of mean absorptions per integration period
        if isempty(cm.absorptions)
            warning('no absorption data');
            return;
        end

        % Show the data, with the gamma from the window.
        axisData.data = mean(cm.absorptions, 3);
        imagesc((axisData.data) .^ gam);
        axis off;

        % Preserve the tick labels in real photons
        colormap(gray);  % Shows a numerical value
        cbar = colorbar;

        % set(get(cbar, 'title'), 'string', 'p per frame', 'rotation', 90);
        photons = str2double(get(cbar, 'TickLabels')) .^ (1 / gam);
        photons = num2str(round(photons));
        set(cbar, 'TickLabels', photons);
        axis image;
        title('Absorptions per integration time');

        set(curAx,'UserData',axisData);  % Put it back

    case 'movieabsorptions'
        % Movie in gray scale

        % Could become cla return
        if isempty(cm.absorptions)
            disp('no absorption data');
            return;
        elseif ismatrix(cm.absorptions)
            disp('no absorption time series')
            return;
        else
            % Additional movie arguments may include the video file name, step,
            % and FrameRate
            ieMovie(cm.absorptions, varargin{:});
        end

    case {'hlineabsorptions', 'vlineabsorptions'}
        % Data are stored in the temporary potting window.
        %
        % cmrPlot(obj,'hline absorptions',[roiPoint])

        data = mean(cm.absorptions, 3);

        if isempty(roi), pt = iePointSelect(app,'Select a point.',1);
        else, pt = roi;
        end
        x = ieClip(pt(1), 1, size(data, 1));
        y = ieClip(pt(2), 1, size(data, 2));
        switch lower(plotType(1))
            case 'h'
                c = [1,size(data,2),y,y];
            case 'v'
                c = [x, x, 1,size(data,1)];
        end

        ieNewGraphWin;
        yStr = 'Absorptions per frame';
        if isequal(plotType(1), 'v')
            plot(data(:, x), 'k-', 'LineWidth', 2);
            uData.y = data(:, x);
            uData.x = 1:length(uData.y);
            grid on;
            xlabel('Vertical position (cones)');
            ylabel(yStr);
            set(curAx, 'userdata', data(:, x));
        else
            plot(data(y, :), 'k-', 'LineWidth', 2);
            uData.y = data(y, :);
            uData.x = 1:length(uData.y);
            grid on;
            xlabel('Horizontal position (cones)');
            ylabel(yStr);
        end

        if isa(app,'coneRectWindow_App')
            ieROIDraw(app,'shape','line','shape data',c);
        end

    case {'hlineabsorptionslms', 'vlineabsorptionslms'}
        % Does not work correctly when in the cone mosaic viewing mode.
        data = mean(cm.absorptions, 3);

        if isempty(roi), pt = iePointSelect(app,'Select a point.',1);
        else, pt = roi;
        end
        x = ieClip(pt(1), 1, size(data, 1));
        y = ieClip(pt(2), 1, size(data, 2));
        switch lower(plotType(1))
            case 'h'
                c = [1,size(data,2),y,y];
            case 'v'
                c = [x, x, 1,size(data,1)];
        end

        if isa(app,'coneRectWindow_App')
            ieROIDraw(app,'shape','line','shape data',c);
        end

        % Get the cone locations in microns
        coneLocs(:,:,1) = reshape(cm.coneLocs(:,1),cm.rows,cm.cols);
        coneLocs(:,:,2) = reshape(cm.coneLocs(:,2),cm.rows,cm.cols);
        coneLocs = coneLocs*1e6;  % In microns

        ieNewGraphWin([], 'tall'); names = 'LMS';
        c = {'ro-', 'go-', 'bo-'};
        yStr = 'Absorptions per frame';
        if isequal(plotType(1), 'v')
            c = {'ro-', 'go-', 'bo-'};
            for ii = 2 : 4 % L, M, S
                % Maybe this should be a coneAbsorptions method?
                subplot(3, 1, ii - 1);
                lst = (cm.pattern(:,x) == ii);
                data = cm.absorptions(lst,x);
                pos = coneLocs(lst,x,2);
                plot(pos, data, c{ii - 1}, 'LineWidth', 2);
                grid on;
                xlabel('Vertical Position (um)');
                ylabel([names(ii - 1) ' ' yStr]);
                uData.pos{ii - 1} = pos;
                uData.data{ii - 1} = data;
            end
        else
            % Horizontal
            for ii = 2:4 % L, M, S
                subplot(3, 1, ii - 1);
                lst = find(cm.pattern(y, :) == ii);
                data = cm.absorptions(y,lst);
                pos = coneLocs(y,lst,1);
                plot(pos, data, c{ii - 1}, 'LineWidth', 2);
                grid on;
                xlabel('Horizontal Position (um)');
                ylabel([names(ii - 1) ' ' yStr]);
                uData.pos{ii - 1} = pos;
                uData.data{ii - 1} = data;
            end
        end

    case 'timeseriesabsorptions'
        % Context menu plot absorption time series.
        data = cm.absorptions;
        if size(data,3) == 1
            warning('No time series.')
            return;
        end

        mx = max(data(:));
        mn = min(data(:));

        if isempty(roi), pt = iePointSelect(app,'Select a point.',1);
        else, pt = roi;
        end
        x = ieClip(pt(1), 1, size(data, 1));
        y = ieClip(pt(2), 1, size(data, 2));
        c = [7 y x];
        if isa(app,'coneRectWindow_App')
            ieROIDraw(app,'shape','circle','shape data',c);
        end

        t = (1:size(data, 3)) * cm.integrationTime * 1e3;

        ieNewGraphWin;
        yStr = 'Absorptions per frame';
        data = squeeze(data(y, x, :));
        plot(t, squeeze(data), 'LineWidth', 2);
        uData.timerseries = t;
        uData.x = t; uData.y = data; uData.pos = [x, y];
        grid on; xlabel('Time (ms)'); ylabel(yStr);
        set(curAx, 'ylim', [mn mx]);

    case 'meancurrent'
        if isempty(cm.current)
            warning('no photocurrent data'); 
            return;
        end

        data = mean(cm.current, 3);

        % Apply gamma. The current is always negative.        
        gam = str2double(app.editGam.Value);
        if max(data(:)) > 0 && ~isequal(gam,1)
            warning('Current is +/-, gamma correction turned off');
            gam = 1;
        end

        % Carry on assuming current is all negative pA.
        % uData = -1*(abs(uData).^gam);
        data = abs(data);
        if ~isequal(hf, 'none'), imagesc(data .^ gam); end
        uData.data = data;        

        % Show the data, with the gamma from the window.
        axisData.data = mean(cm.current, 3);
        if isequal(gam,1), imagesc(axisData.data);
        else,              imagesc((axisData.data) .^ gam);
        end
        axis off;

        % Preserve the tick labels in real photons
        colormap(flipud(gray));  % Shows a numerical value
        cbar = colorbar;
        current = -1 * ...
            (abs(str2double(get(cbar, 'TickLabels')) .^ (1 / gam)));
        current = num2str(round(current));
        set(cbar, 'TickLabels', current);
        axis image;
        title('Photocurrent (pA)');

        set(curAx,'UserData',axisData);  % Put it back

    case {'hlinecurrent', 'vlinecurrent'}
        data = mean(cm.current, 3);

        pt = iePointSelect(app,'Select a point.',1);
        x = ieClip(pt(1), 1, size(data, 1));
        y = ieClip(pt(2), 1, size(data, 2));
        switch lower(plotType(1))
            case 'h'
                c = [1,size(data,2),y,y];
            case 'v'
                c = [x, x, 1,size(data,1)];
        end
        if isa(app,'coneRectWindow_App')
            ieROIDraw(app,'shape','line','shape data',c);
        end

        ieNewGraphWin;
        yStr = 'Absorptions per frame';
        if isequal(plotType(1), 'v')
            plot(data(:, x), 'LineWidth', 2);
            grid on;
            xlabel('Vertical position (cones)');
            ylabel(yStr);
            uData.y = data(:, x); uData.x = 1:length(uData.y);
        else
            plot(data(y, :), 'LineWidth', 2);
            grid on;
            xlabel('Horizontal position (cones)');
            ylabel(yStr);
            uData.y = data(y, :);
            uData.x = 1:length(uData.y);
        end

    case {'hlinecurrentlms', 'vlinecurrentlms'}
        data = mean(cm.current, 3);

        pt = iePointSelect(app,'Select a point.',1);
        x = ieClip(pt(1), 1, size(data, 1));
        y = ieClip(pt(2), 1, size(data, 2));
        switch lower(plotType(1))
            case 'h'
                c = [1,size(data,2),y,y];
            case 'v'
                c = [x, x, 1,size(data,1)];
        end
        if isa(app,'coneRectWindow_App')
            ieROIDraw(app,'shape','line','shape data',c);
        end

        ieNewGraphWin([], 'tall');
        names = 'LMS';
        c = {'ro-', 'go-', 'bo-'};
        yStr = 'Photocurrent (pA)';
        if isequal(plotType(1), 'v')
            c = {'ro-', 'go-', 'bo-'};
            for ii = 2:4 % L, M, S
                subplot(3, 1, ii - 1);
                pos = find(cm.pattern(:, x) == ii);
                plot(pos, data(pos, x), c{ii-1}, 'LineWidth', 2);
                grid on;
                uData.pos{ii - 1} = pos;
                uData.data{ii - 1} = data(pos, x);
                xlabel('Vertical position (cones');
                ylabel([names(ii-1) ' ' yStr]);
                set(curAx, 'xlim', [1 size(data, 1)]);
            end
        else
            for ii = 2:4 % L, M, S
                subplot(3, 1, ii - 1);
                pos = find(cm.pattern(y, :) == ii);
                plot(pos, data(y, pos), c{ii - 1}, 'LineWidth', 2);
                grid on;
                uData.pos{ii - 1} = pos;
                uData.data{ii - 1} = data(y, pos);
                xlabel('Horizontal position (cones');
                ylabel([names(ii - 1) ' ' yStr]);
                set(curAx, 'xlim', [1 size(data, 2)]);
            end
        end
    case 'timeseriescurrent'
        data = cm.current;
        mx = max(data(:));
        mn = min(data(:));

        pt = iePointSelect(app,'Select a point.',1);
        x = ieClip(pt(1), 1, size(data, 1));
        y = ieClip(pt(2), 1, size(data, 2));
        c = [3 x y];
        if isa(app,'coneRectWindow_App')
            ieROIDraw(app,'shape','circle','shape data',c);
        end

        t = (1:size(data, 3)) * cm.integrationTime * 1e3;

        ieNewGraphWin;
        yStr = 'Absorptions per frame';
        plot(t, squeeze(data(y, x, :)), 'LineWidth', 2);
        uData.x = t;
        uData.y = squeeze(data(y, x, :));
        uData.pos = [x, y];
        grid on;
        xlabel('Time (ms)');
        ylabel(yStr);
        set(curAx, 'ylim', [mn mx]);

    case 'impulseresponse'
        % Current impulse response at cone mosaic temporal sampling rate
        %
        % Not called from within the window, so it should probably be
        % moved out of here.


        % The outersegment cone temporal impulse response functions are
        % always represented at a high sampling rate (0.1 ms).

        lmsFilters = cm.os.linearFilters(obj);

        % Old code if ~isempty(cm.absorptions)
        % This doesn't make sense to me (BW).  Guessing this
        % if/else is older code that should be removed.
        %{
            absorptionsInXWFormat = RGB2XWFormat(cm.absorptions);
            lmsFilters = cm.os.linearFilters('absorptionsInXWFormat', ...
                absorptionsInXWFormat);
        %}



        %% Interpolate stored lmsFilters to the time base of absorptions
        osTimeAxis = cm.os.timeAxis;
        coneTimeAxis = cm.interpFilterTimeAxis;

        % Interpolate down to the cone mosaic sampling rate Interpolation
        % assumes that we are accounting for the time sample bin width
        % elsewhere. Also, we extrapolate the filters with zeros to make
        % sure that they extend all the way through the absorption time
        % axis. See the notes in s_matlabConv2.m for an explanation of why.
        interpFilters = interp1(osTimeAxis(:), lmsFilters, ...
            coneTimeAxis(:), 'linear', 0);

        plot(coneTimeAxis, interpFilters(:, 1), 'r-o', ...
            coneTimeAxis, interpFilters(:, 2), 'g-o', ...
            coneTimeAxis, interpFilters(:, 3), 'b-o');
        xlabel('Time (sec)');
        ylabel('Current (pA)');
        grid on;
        l = {'L cone', 'M cone', 'S cone'};
        legend(l);
        title('Impulse response (cone temporal sampling)');

    case 'moviecurrent'
        % Current movie in gray scale
        if isempty(cm.current)
            if isempty(p.Results.hf), close(hf); end
            error('no current data');
        end

        % Additional arguments may be the video file name, step, and
        % FrameRate
        ieMovie(cm.current, varargin{:});

    case 'conefundamentals'
        % The cone absorptance without macular pigment or lens
        uData = cm.pigment.absorptance;
        if ~isequal(hf, 'none')
            plot(cm.wave, uData, 'LineWidth', 2);
            grid on;
            xlabel('Wavelength (nm)');
            ylabel('Cone quanta absorptance');
        end

    case 'maculartransmittance'
        uData = cm.macular.transmittance;
        if ~isequal(hf, 'none')
            plot(cm.wave, uData, 'LineWidth', 2);
            grid on;
            xlabel('Wavelength (nm)');
            ylabel('Macular transmittance');
        end

    case 'macularabsorptance'
        uData = cm.macular.absorptance;
        if ~isequal(hf, 'none')
            plot(cm.wave, uData, 'LineWidth', 2);
            grid on;
            xlabel('Wavelength (nm)');
            ylabel('Macular absorptance');
        end

    case 'macularabsorbance'
        uData = cm.macular.unitDensity;
        if ~isequal(hf, 'none')
            plot(cm.wave, uData, 'LineWidth', 2);
            grid on;
            xlabel('Wavelength (nm)');
            ylabel('Macular absorbance');
        end

    case 'conespectralqe'
        % Quantum efficiency of macular pigment and cone photopigments
        uData = cm.qe;
        if ~isequal(hf, 'none')
            plot(cm.wave, cm.qe, 'LineWidth', 2);
            grid on;
            xlabel('Wavelength (nm)');
            ylabel('Cone quanta efficiency');
        end

    case 'eyespectralqe'
        % Includes lens, macular pigment, and cone photopigment properties
        if isempty(oi), error('oi required for spectral qe'); end
        lensTransmittance = oiGet(oi, 'lens transmittance', ...
            'wave', cm.wave);
        uData = bsxfun(@times, lensTransmittance, cm.qe);

        if ~isequal(hf, 'none')
            plot(cm.wave, uData, 'LineWidth', 2);
            grid on;
            xlabel('Wavelength (nm)');
            ylabel('Eye quanta efficiency');
        end

        % case 'currenttimeseries'
        %     % Photocurrent time series of selected points.
        %     % Need a way to choose which points!
        %     if isempty(cm.current)
        %         if isempty(p.Results.hf), close(hf); end
        %         error('no photocurrent data');
        %     end
        %     uData = plotCurrentTimeseries(obj, varargin{:});

    case {'empath', 'eyemovementpath'}
        % coneRectPlot(app,'eye movement path');

        plot(cm.emPositions(:, 1), cm.emPositions(:, 2),'ko:');
        xLim = [min(cm.emPositions(:,1)),max(cm.emPositions(:,1))];
        yLim = [min(cm.emPositions(:,2)),max(cm.emPositions(:,2))];
        if xLim(1) > -1, xLim(1) = -3; end; if xLim(2) < 1,  xLim(2) = 3; end
        if yLim(1) > -1, yLim(1) = -3; end; if yLim(2) < 1,  yLim(2) = 3; end
        grid on;
        xlabel('Horizontal position (cones)');
        ylabel('Vertical position (cones)');
        set(curAx,'xlim',xLim,'ylim',yLim);
        title(sprintf('Eye movement path (%.1f ms steps)',cm.integrationTime*1e3));

        % RGB movies on cone mosaic. These are not currently implemented,
        % but exist here in draft form. See routine coneImageActivity below
        % as well. Could be resurrected some day.
        %
        % case 'absorptions'
        %     % Movie of the absorptions on the cone mosaic
        %     if isempty(cm.absorptions)
        %         % Could be come cla; return
        %         error('no absorption data');
        %     end
        %     uData = coneImageActivity(obj, hf, varargin{:});
        %
        % case {'current', 'photocurrent'}
        %     % Photo current movie on colored cone mosaic
        %     if isempty(cm.current)
        %         if isempty(p.Results.hf), close(hf); end
        %         error('no photocurrent data');
        %     end
        %     uData = coneImageActivity(obj, hf, 'dataType', ...
        %         'photocurrent', varargin{:});

    otherwise
        error('unsupported plot type');
end

% Put back the modified user data
if exist('uData','var'), set(curAx, 'userdata', uData); end

end

%{
function mov = coneImageActivity(cMosaic, hf, varargin)
% Make a movie or a single image of cone absorptions on a colored mosaic
%
% Syntax:
%   mov = coneImageActivity(coneMosaicRect, hf)
%
% Description:
%    This function Would be used in commented out cases above if they are
%    resurrected.
%
%    Make a movie or a single image of cone absorptions on a colored mosaic
%
% Inputs:
%    cones - coneMosaicRect class object=
%    hf    - figure handle, 'none', a boolean, or a struct.
%        Struct  - If a struct, hf contains the movie parameters. The movie
%                  is saved based on the name field in these parameters.
%                  The parameters are .vname (video file name) and
%                  .FrameRate (video frame rate).
%        Boolean - If dFLag is a boolean, the value indicates whether to
%                  show the movie (true) or not (false).
%
% Outputs:
%    mov   - The created movie.
%
%  Optional key/value pairs:
%    **Fill out**
%
% Notes:
%    * [Note: DHB - SAY WHAT FORMAT THE RETURNED MOVIE IS IN.]
%    * [Note: DHB - THERE ARE A BUNCH OF KEY/VALUE PAIRS DEFINED BUT NOT
%      DOCUMENTED HERE.]
%

% History:
%    xx/xx/16  HJ/BW  ISETBIO Team, 2016
%    02/19/18  jnm    Formatting

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('cMosaic', @(x) (isa(x, 'coneMosaicRect')));
p.addRequired('hf', @(x) (ischar(x) || ishandle(x) || isstruct(x)));
p.addParameter('step', [], @isnumeric);
p.addParameter('showBar', true, @islogical);
p.addParameter('gamma', 0.3, @isnumeric);
p.addParameter('dataType', 'absorptions', @ischar);

% set parameters
p.parse(cMosaic, hf, varargin{:});
step = p.Results.step;
showBar = p.Results.showBar;
gamma = p.Results.gamma;

switch ieParamFormat(p.Results.dataType)
    case {'absorptions', 'isomerizations', 'photons'}
        data = cMosaic.absorptions;
    case {'current', 'photocurrent'}
        data = cMosaic.current;
    otherwise
        error('Unknown data type');
end

nframes = size(data, 3);
if isempty(step), step = max(1, nframes / 100); end

% create cone mosaic image
[~, mosaicData] = cMosaic.plot('cone mosaic', 'hf', 'none');
coneMosaicImage = mosaicData.mosaicImage;
delta = mosaicData.delta;

% bulid frame by frame
if showBar, wbar = waitbar(0, 'creating cone movie'); end
mov = zeros([size(coneMosaicImage), length(1:step:nframes)]);
g = fspecial('gaussian', 6, 2);

for ii = 1:step:nframes
    if showBar, waitbar(ii / nframes, wbar); end

    % Expand the data to the size of the coneMosaicImage
    d = cMosaic.absorptions(:, :, ii);
    fgrid = full(ffndgrid(cMosaic.coneLocs * 1e6, d(:), delta));

    % Scale the cone mosaic image by the adapted data
    index = (ii - 1) / step + 1;
    for jj = 1:3
        mov(:, :, jj, index) = convolvecirc(coneMosaicImage(:, :, jj) ...
            .* fgrid, g);
    end
end

% Scale and gamma correct mov
mov = ieScale(mov) .^ gamma;

if showBar, close(wbar); end

% show the movie, or write to file
if isstruct(hf)
    % When dFlag is a struct, show the move and save it in a file
    vObj = VideoWriter(hf.vname);
    vObj.FrameRate = hf.FrameRate;
    open(vObj);
    for ii = 1:size(mov, 4)
        image(mov(:, :, :, ii));
        drawnow;
        F = getframe;
        writeVideo(vObj, F);
    end
    close(vObj);
elseif ~isequal(hf, 'none') && ishandle(hf)
    % If it is a figure handle, show it in that figure
    for ii = 1:size(mov, 4), imshow(mov(:, :, :, ii)); drawnow; end
end

end

function uData = plotCurrentTimeseries(obj, varargin)
% Pull out the time series of the photo current and plot it
%
% Syntax:
%   uData = plotCurrentTimeseries(obj, [varargin])
%
% Description:
%    Pull out the time series of the photo current and plot it
%
% Inputs:
%    obj   - The cone mosaic object
%
% Outputs:
%    uData - The user data
%
% Optional key/value pairs:
%    Need to fill out!
%

% Temporal samples
dt = cm.os.timeStep;

% The current
outputSignalTemp = cm.current;
sz = size(outputSignalTemp);

% Reshape for plotting
outputSignal = reshape(outputSignalTemp, sz(1) * sz(2), sz(3));

% Plot a random subset ... must handle more explicitly
uData.time = (0:size(outputSignal, 2) - 1) * dt;
uData.current = ...
    outputSignal(1 + floor((sz(1) * sz(2) / 100) * rand(200, 1)), :);
plot(uData.time, uData.current);

% title('Output current');
xlabel('Time (sec)');
ylabel('pA');

end
%}

%----------------------
function coneRectPlotHelp(validPlotsPrint)

fprintf('\nconeRectPlot plot types\n--------------\n');
for ii = 2:length(validPlotsPrint), fprintf('\t%s\n', validPlotsPrint{ii}); end

end