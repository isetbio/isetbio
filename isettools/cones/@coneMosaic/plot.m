function [uData, hf] = plot(obj, plotType, varargin)
%PLOT  Plot function for @conemmsaic base class
%
%   [uData, hf] = plot(obj,plotType, varargin)
%
% There is a specialized plot method for the coneMosaicHex class that calls
% this function.
% 
% When the plot type string begins with 'os ' or 'outersegment ' we pass
% the arguments along to os.plot().  For example,
%
%      cMosaic.plot('os impulse response')
%
% graphs the outer segment impulse response on its own time axis.
% [DHB NOTE: RETURN DATA ARE BOTH SET TO EMPTY IN THIS CASE, NOT GOOD.
% CAN THIS BE EASILY FIXED?]
%
% Required Inputs:
%   plotType - string, type of plot
%
% Optional key/value pairs:
%  'hf' - figure handle or control structure, the meaning of value is
%     [] - create plot in new figure
%    'none'        - don't plot
%     isgraphics figure or axes handle - Use that figure or axes
%  'oi' - an optical image that is used to get lens transmittance for QE
%         plots 
%  'x' - x value (col) for vline plots 
%  'y' - y value (row) for hline plots
%   
% Outputs:
%   hf    - figure or axes handle
%   uData - Computed data
%
% @conemosaic.plot('help') lists the valid plot types
%  N.B. The '*' means user must select a point on an image to specify the
%  location of the point or line to be plotted).
%
%   Cone mosaic            - Color image of the cone arrangement
%
% 	Mean absorptions       - Image of the mean absorptions
% 	Movie absorptions      - Movie of the absorptions on cone mosaic
% 	Mean current           - Image of the mean current
% 	Movie current          - Gray scale movie of current
%
% 	hline absorptions*     - Graph of a horizontal line of absoprtions
% 	hline current*         - 
% 	hline absorptions lms* - Three panel graph of LMS absorptions
% 	hline current lms*     - Three panel graph of LMS current
% 	vline absorptions*     - Vertical line
% 	vline current*
% 	vline absorptions lms*
% 	vline current lms*
% 	time series absorptions* - Cone absorptions Graph of a point
% 	time series current*     - Cone photocurrent Graph of a point
%
% 	Impulse response       - Cone current impulse response
%
% 	Cone fundamentals      - Cone pigment without macular or lens
% 	Cone spectral QE       - Cone pigment and macular
% 	Eye spectral QE        - Cone pigment with macular and lens
% 	Macular transmittance  - Graph
% 	Macular absorptance    - Graph
% 	Macular absorbance     - Graph
% 	Eye movement path      - eye movement path
%     
%
%  Examples:
%   coneMosaic.plot('impulse response')
%   coneMosaic.plot('cone mosaic')
%
% HJ/BW, ISETBIO TEAM, 2016

% Programming TODO
% Think about coneImageActivity function at end.

%% Check plot type string if we send this off to the os plot routine
%
% (Might Find a cleaner way to check and send to os.plot.
%  Maybe create a parse argument string as in ISET.)
if (length(plotType) > 3 && strcmp(plotType(1:3),'os '))
    obj.os.plot(plotType(4:end),'cmosaic',obj,varargin{:});
    return;
elseif (length(plotType) > 13 && strcmp(plotType(1:13),'outersegment '))
    obj.os.plot(plotType(14:end),'cmosaic',obj,varargin{:});
    return;
end

%% It is a cone mosaic plot, not an os plot.

% Add to this list when you put in a new type of plot.  The validation
% squeezes out the spaces and makes lower case.
validPlots = {'help',...
    'Cone mosaic',...
    'Mean absorptions','Movie absorptions', ...
    'Mean current','Movie current',...
    'hline absorptions','hline current','hline absorptions lms','hline current lms',...
    'vline absorptions','vline current','vline absorptions lms', 'vline current lms' ...
    'Impulse response','time series absorptions','time series current'...
    'Cone fundamentals', 'Cone spectral QE', 'Eye spectral QE',...
    'Macular transmittance','Macular absorptance', 'Macular absorbance', ...
    'Eye movement path'};

% If the person just wants help
if isequal(plotType,'help')
    fprintf('\nKnown %s plot types\n--------------\n',class(obj));
    for ii=2:length(validPlots)
        fprintf('\t%s\n',validPlots{ii});
    end
    return;
else
    % Makes less readable, but better for programming. I
    for ii=1:length(validPlots)
        validPlots{ii} = ieParamFormat(validPlots{ii});
    end
end


%% Onward
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('obj');

% Squeeze spaces and force lower case on validPlots.  Then check plotType
validPlots = cellfun(@(x)(ieParamFormat(x)),validPlots,'UniformOutput',false);
p.addRequired('pType', @(x)any(validatestring(ieParamFormat(x),validPlots)));

p.addParameter('hf', []);                     % Figure handle
p.addParameter('oi',[],@isstruct);            % Used for spectral qe
p.addParameter('x',[],@isscalar);             % x axis value
p.addParameter('y',[],@isscalar);             % y axis value

p.parse(obj,plotType, varargin{:});
hf = p.Results.hf;
oi = p.Results.oi;                         % Used in plotGraphs routine

% Initialize return structure
uData = [];

%% Initialize where we'll plot
if isempty(hf), hf = vcNewGraphWin;
elseif isgraphics(hf, 'figure'), figure(hf);
elseif isgraphics(hf, 'axes'), axes(hf);
end

%% Set color order so that LMS plots as RGB
if ~isequal(hf, 'none')
    co = get(gca, 'ColorOrder');
    if isgraphics(hf,'axes')
        set(get(hf,'parent'),'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :))
    else
        set(hf, 'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :));
    end
end

%% Switch on passed plot type
switch ieParamFormat(plotType)
    
    case 'conemosaic'
        % Default for cone size
        support = [4,4]; spread = 2; maxCones = 5e4;
        
        % Speed things up
        nCones = size(obj.coneLocs,1);
        locs = obj.coneLocs; pattern = obj.pattern(:);
        if  nCones > maxCones
            disp('Displaying subsampled (50K) version')
            lst = randi(nCones,[maxCones,1]);
            lst = unique(lst);
            locs = locs(lst,:); pattern = pattern(lst,:);
            
            % Brighten up in this case
            support = round([nCones/maxCones,nCones/maxCones]);
            spread = 2*support(1);
        end
        
        [uData.support, uData.spread, uData.delta, uData.mosaicImage] = ...
            conePlot(locs * 1e6, pattern, support, spread);
        imagesc(uData.mosaicImage); axis off; axis image;
        
    case 'meanabsorptions'
        % Image of mean absorptions per integration period
        %
        if isempty(obj.absorptions), error('no absorption data'); end
        
        % Show the data, with the gamma from the window.
        uData.data = mean(obj.absorptions,3);
        gdata      = guidata(obj.hdl);
        gam        = str2double(get(gdata.editGam,'string'));
        imagesc((uData.data).^gam); axis off;
        
        % Preserve the tick labels in real photons
        colormap(gray);  % Shows a numerical value
        cbar = colorbar;
        
        % set(get(cbar,'title'),'string','p per frame','rotation',90);
        photons = str2double(get(cbar,'TickLabels')).^(1/gam);
        photons = num2str(round(photons)); set(cbar,'TickLabels',photons);
        axis image; title('Absorptions per integration time');
        
    case 'movieabsorptions'
        % Movie in gray scale
        
        % Could become cla return
        if isempty(obj.absorptions), error('no absorption data'); end
        
        % Additional arguments may be the video file name, step, and
        % FrameRate
        uData = ieMovie(obj.absorptions,varargin{:});
        
    case {'hlineabsorptions','vlineabsorptions'}
        data = mean(obj.absorptions,3);
        
        % The plots below are with respect to a point.
        % Get the point
        [x, y] = ginput(1); % Rounded and clipped to the data
        x = ieClip(round(x), 1, size(data, 2));
        y = ieClip(round(y), 1, size(data, 1));
        
        % Draw a circle around the selected point.
        viscircles([x,y],0.7);
        vcNewGraphWin;
        yStr = 'Absorptions per frame';
        if isequal(plotType(1),'v')
            plot(data(:, x), 'LineWidth', 2);
            uData.y = data(:,x); uData.x = 1:length(uData.y);
            grid on; xlabel('Vertical position (cones)'); ylabel(yStr);
            set(gca,'userdata',data(:,x));
        else
            plot(data(y, :), 'LineWidth', 2);
            uData.y = data(y, :); uData.x = 1:length(uData.y);
            grid on; xlabel('Horizontal position (cones)'); ylabel(yStr);
        end

    case {'hlineabsorptionslms','vlineabsorptionslms'}
        data = mean(obj.absorptions,3);
        
        % The plots below are with respect to a point.
        % Get the point
        [x, y] = ginput(1); % Rounded and clipped to the data
        x = ieClip(round(x), 1, size(data, 2));
        y = ieClip(round(y), 1, size(data, 1));        %
        viscircles([x,y],0.7);

        vcNewGraphWin([],'tall'); names = 'LMS';
        c = {'ro-','go-','bo-'};
        yStr = 'Absorptions per frame';
        if isequal(plotType(1),'v')
            c = {'ro-','go-','bo-'};
            for ii = 2 : 4 % L, M, S
                subplot(3, 1, ii-1);
                pos = find(obj.pattern(:, x) == ii);
                plot(pos, data(pos, x), c{ii-1}, 'LineWidth', 2); grid on;
                uData.pos{ii-1} = pos; uData.data{ii-1}=data(pos,x);
                xlabel('Vertical Position (cones');
                ylabel([names(ii-1) ' ' yStr]);
                set(gca,'xlim',[1 size(data,1)]);
            end
        else
            for ii = 2 : 4 % L, M, S
                subplot(3, 1, ii-1);
                pos = find(obj.pattern(y, :) == ii);
                plot(pos, data(y, pos), c{ii-1}, 'LineWidth', 2); grid on;
                uData.pos{ii-1} = pos; uData.data{ii-1}=data(y,pos);
                xlabel('Horizontal Position (cones');
                ylabel([names(ii-1) ' ' yStr]);
                set(gca,'xlim',[1 size(data,2)]);
            end
        end

    case 'timeseriesabsorptions'
        % Context menu plot absorption time series.
        data = obj.absorptions;
        mx = max(data(:));
        mn = min(data(:));
                
        [x, y] = ginput(1); % Rounded and clipped to the data
        x = ieClip(round(x), 1, size(data, 2));
        y = ieClip(round(y), 1, size(data, 1));        %
        viscircles([x,y],0.7);
        
        t = (1:size(data, 3)) * obj.integrationTime * 1e3;

        vcNewGraphWin;         
        yStr = 'Absorptions per frame';
        data = squeeze(data(y, x, :));
        plot(t, squeeze(data), 'LineWidth', 2);
        uData.timerseries = t;
        uData.x = t; uData.y = data; uData.pos = [x,y];
        grid on; xlabel('Time (ms)'); ylabel(yStr);
        set(gca,'ylim',[mn mx]);
        
    case 'meancurrent'
        if isempty(obj.current), error('no photocurrent data'); end
        data = mean(obj.current, 3);
        
        % Apply gamma.  The current is always negative.
        gdata = guidata(obj.hdl);
        gam = str2double(get(gdata.editGam,'string'));
        if max(data(:)) > 0
            warning('Gamma correction in display is not correct');
        end
        
        % Carry on assuming current is negative pA.
        % uData = -1*(abs(uData).^gam);
        data = abs(data);
        if ~isequal(hf, 'none'), imagesc(data.^gam); end
        uData.data = data;
        
        axis off; colormap(flipud(gray));  % Shows a numerical value
        cbar = colorbar;
        current = -1*(abs(str2double(get(cbar,'TickLabels')).^(1/gam)));
        current = num2str(round(current)); set(cbar,'TickLabels',current);
        axis image;
        title('Photocurrent (pA)');
        
    case {'hlinecurrent','vlinecurrent'}
        data = mean(obj.current,3);
        
        % The plots below are with respect to a point.
        % Get the point
        [x, y] = ginput(1); % Rounded and clipped to the data
        x = ieClip(round(x), 1, size(data, 2));
        y = ieClip(round(y), 1, size(data, 1));
        
        % Draw a circle around the selected point.
        viscircles([x,y],0.7);
        vcNewGraphWin;
        yStr = 'Absorptions per frame';
        if isequal(plotType(1),'v')
            plot(data(:, x), 'LineWidth', 2);
            grid on; xlabel('Vertical position (cones)'); ylabel(yStr);
            uData.y = data(:,x); uData.x = 1:length(uData.y); 
        else
            plot(data(y, :), 'LineWidth', 2);
            grid on; xlabel('Horizontal position (cones)'); ylabel(yStr);
            uData.y = data(y,:);
            uData.x = 1:length(uData.y);
        end
        
    case {'hlinecurrentlms','vlinecurrentlms'}
        
        data = mean(obj.current,3);
        
        % The plots below are with respect to a point.
        % Get the point
        [x, y] = ginput(1); % Rounded and clipped to the data
        x = ieClip(round(x), 1, size(data, 2));
        y = ieClip(round(y), 1, size(data, 1));        %
        viscircles([x,y],0.7);

        vcNewGraphWin([],'tall'); names = 'LMS';
        c = {'ro-','go-','bo-'};
        yStr = 'Photocurrent (pA)';
        if isequal(plotType(1),'v')
            c = {'ro-','go-','bo-'};
            for ii = 2 : 4 % L, M, S
                subplot(3, 1, ii-1);
                pos = find(obj.pattern(:, x) == ii);
                plot(pos, data(pos, x), c{ii-1}, 'LineWidth', 2); grid on;
                uData.pos{ii-1} = pos; uData.data{ii-1}=data(pos,x);
                xlabel('Vertical Position (cones');
                ylabel([names(ii-1) ' ' yStr]);
                set(gca,'xlim',[1 size(data,1)]);
            end
        else
            for ii = 2 : 4 % L, M, S
                subplot(3, 1, ii-1);
                pos = find(obj.pattern(y, :) == ii);
                plot(pos, data(y, pos), c{ii-1}, 'LineWidth', 2); grid on;
                uData.pos{ii-1} = pos; uData.data{ii-1}=data(y,pos);
                xlabel('Horizontal Position (cones');
                ylabel([names(ii-1) ' ' yStr]);
                set(gca,'xlim',[1 size(data,2)]);
            end
        end
    case 'timeseriescurrent'
        data = obj.current;
        mx = max(data(:));
        mn = min(data(:));
        
        [x, y] = ginput(1); % Rounded and clipped to the data
        x = ieClip(round(x), 1, size(data, 2));
        y = ieClip(round(y), 1, size(data, 1));        %
        viscircles([x,y],0.7);
        
        t = (1:size(data, 3)) * obj.integrationTime * 1e3;
        
        vcNewGraphWin;
        yStr = 'Absorptions per frame';
        plot(t, squeeze(data(y, x, :)), 'LineWidth', 2);
        uData.x = t; uData.y = squeeze(data(y, x, :));
        uData.pos = [x, y];
        grid on; xlabel('Time (ms)'); ylabel(yStr);
        set(gca,'ylim',[mn mx]);
        
    case 'impulseresponse'
        % The current impulse response at the cone mosaic temporal sampling
        % rate
        
        % The outersegment cone temporal impulse response functions are
        % always represented at a high sampling rate (0.1 ms).
        if isempty(obj.absorptions)
            lmsFilters = obj.os.linearFilters(obj);
        else
            absorptionsInXWFormat = RGB2XWFormat(obj.absorptions);
            lmsFilters = obj.os.linearFilters('absorptionsInXWFormat', absorptionsInXWFormat);
        end
        
        %% Interpolate the stored lmsFilters to the time base of the absorptions
        osTimeAxis    = obj.os.timeAxis;
        coneTimeAxis  = obj.interpFilterTimeAxis;
        
        % Interpolate down to the cone mosaic sampling rate
        % Interpolation assumes that we are accounting for the time sample bin
        % width elsewhere.  Also, we extrapolate the filters with zeros to make
        % sure that they extend all the way through the absorption time axis.
        % See the notes in s_matlabConv2.m for an explanation of why.
        interpFilters = interp1(osTimeAxis(:),lmsFilters,coneTimeAxis(:),'linear',0);
        
        vcNewGraphWin;
        plot(coneTimeAxis,interpFilters(:,1),'r-', ...
            coneTimeAxis,interpFilters(:,2),'g-', ...
            coneTimeAxis,interpFilters(:,3),'b-o');
        xlabel('Time (sec)'); ylabel('Current (pA)');
        grid on;
        l = {'L cone','M cone','S cone'}; legend(l);
        title('Impulse response (cone temporal sampling)');
        
    case 'moviecurrent'
        % Current movie in gray scale
        if isempty(obj.current)
            if isempty(p.Results.hf), close(hf); end
            error('no current data');
        end
        
        % Additional arguments may be the video file name, step, and
        % FrameRate
        uData = ieMovie(obj.current,varargin{:});
        
    case 'conefundamentals'
        % The cone absorptance without macular pigment or lens
        uData = obj.pigment.absorptance;
        if ~isequal(hf, 'none')
            plot(obj.wave, uData, 'LineWidth', 2); grid on;
            xlabel('Wavelength (nm)'); ylabel('Cone quanta absorptance');
        end
        
    case 'maculartransmittance'
        uData = obj.macular.transmittance;
        if ~isequal(hf, 'none')
            plot(obj.wave, uData, 'LineWidth', 2); grid on;
            xlabel('Wavelength (nm)'); ylabel('Macular transmittance');
        end
        
    case 'macularabsorptance'
        uData = obj.macular.absorptance;
        if ~isequal(hf, 'none')
            plot(obj.wave, uData, 'LineWidth', 2); grid on;
            xlabel('Wavelength (nm)'); ylabel('Macular absorptance');
        end
        
    case 'macularabsorbance'
        uData = obj.macular.unitDensity;
        if ~isequal(hf, 'none')
            plot(obj.wave, uData, 'LineWidth', 2); grid on;
            xlabel('Wavelength (nm)'); ylabel('Macular absorbance');
        end
        
    case 'conespectralqe'
        % Quantum efficiency of macular pigment and cone photopigments
        uData = obj.qe;
        if ~isequal(hf, 'none')
            plot(obj.wave, obj.qe, 'LineWidth', 2); grid on;
            xlabel('Wavelength (nm)'); ylabel('Cone quanta efficiency');
        end
        
    case 'eyespectralqe'
        % Includes lens, macular pigment, and cone photopigment properties
        if isempty(oi), error('oi required for spectral qe'); end
        lensTransmittance = oiGet(oi,'lens transmittance','wave',obj.wave);
        uData = bsxfun(@times, lensTransmittance, obj.qe);
        
        if ~isequal(hf, 'none')
            plot(obj.wave, uData, 'LineWidth', 2); grid on;
            xlabel('Wavelength (nm)'); ylabel('Eye quanta efficiency');
        end
        
        %     case 'currenttimeseries'
        %         % Photocurrent time series of selected points.
        %         % Need a way to choose which points!
        %         if isempty(obj.current)
        %             if isempty(p.Results.hf), close(hf); end
        %             error('no photocurrent data');
        %         end
        %         uData = plotCurrentTimeseries(obj,varargin{:});
        
    case {'empath', 'eyemovementpath'}
        plot(obj.emPositions(:, 1), obj.emPositions(:, 2));
        grid on; xlabel('Horizontal position (cones)');
        ylabel('Vertical position (cones)');
        
        % RGB movies on cone mosaic. These are not currently implemented, but exist
        % here in draft form.  See routine coneImageActivity below as well.
        % Could be resurrected some day.
        %
        %     case 'absorptions'
        %         % Movie of the absorptions on the cone mosaic
        %         if isempty(obj.absorptions)
        %             % Could be come cla; return
        %             error('no absorption data');
        %         end
        %         uData = coneImageActivity(obj, hf, varargin{:});
        %
        %     case {'current', 'photocurrent'}
        %         % Photo current movie on colored cone mosaic
        %         if isempty(obj.current)
        %             if isempty(p.Results.hf), close(hf); end
        %             error('no photocurrent data');
        %         end
        %         uData = coneImageActivity(obj, hf, 'dataType', ...
        %             'photocurrent', varargin{:});
        
        
    otherwise
        error('unsupported plot type');
end

set(gca,'userdata',uData);

end

function mov = coneImageActivity(cMosaic, hf, varargin)
%CONEIMAGEACTIVITY  Make a movie or a single image of cone absorptions on a colored mosaic
%  mov = coneImageActivity(coneMosaic,hf)
%
% This function Would be used in commented out cases above if they are resurrected.
%
%  Inputs:
%  cones:  coneMosaic class object=
%  hf:     figure handle or 'none' or a struct
%          If a struct, hf contains the movie parameters.  The movie is
%          saved based on the name field in these parameters. The parameters
%          are .vname (video file name) and .FrameRate (video frame rate).
%          If dFLag is a boolean, the value indicates whether to show the
%          movie (true) or not (false).
%
%  Outputs:
%  mov     [DHB NOTE: SAY WHAT FORMAT THE RETURNED MOVIE IS IN.]
%
%  Optional key/value pairs:
%    [DHB NOTE: THERE ARE A BUNCH OF KEY/VALUE PAIRS DEFINED BUT NOT DOCUMENTED HERE.]

% HJ/BW, ISETBIO Team, 2016

%% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('cMosaic',@(x) (isa(x,'coneMosaic')));
p.addRequired('hf',@(x) (ischar(x) || ishandle(x) || isstruct(x)));
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
if isempty(step), step = max(1, nframes/100); end

% create cone mosaic image
[~, mosaicData] = cMosaic.plot('cone mosaic', 'hf', 'none');
coneMosaicImage = mosaicData.mosaicImage;
delta = mosaicData.delta;

% bulid frame by frame
if showBar, wbar = waitbar(0,'creating cone movie'); end
mov = zeros([size(coneMosaicImage), length(1:step:nframes)]);
g = fspecial('gaussian', 6, 2);

for ii = 1 : step : nframes
    if showBar, waitbar(ii/nframes,wbar); end
    
    % Expand the data to the size of the coneMosaicImage
    d = cMosaic.absorptions(:, :, ii);
    fgrid = full(ffndgrid(cMosaic.coneLocs * 1e6, d(:), delta));
    
    % Scale the cone mosaic image by the adapted data
    index = (ii - 1)/step + 1;
    for jj = 1 : 3
        mov(:,:,jj,index) = convolvecirc(coneMosaicImage(:,:,jj).*fgrid,g);
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
    for ii = 1:size(mov,4)
        image(mov(:,:,:,ii)); drawnow;
        F = getframe; writeVideo(vObj,F);
    end
    close(vObj);
elseif ~isequal(hf, 'none') && ishandle(hf)
    % If it is a figure handle, show it in that figure
    for ii=1:size(mov,4), imshow(mov(:,:,:,ii)); drawnow; end
end

end


function uData = plotCurrentTimeseries(obj,varargin)
%PLOTCURRENTTIMESEREIS  Pull out the time series of the photo current and plot it

% Temporal samples
dt = obj.os.timeStep;

% The current
outputSignalTemp = obj.current;
sz = size(outputSignalTemp);

% Reshape for plotting
outputSignal = reshape(outputSignalTemp,sz(1)*sz(2),sz(3));

% Plot a random subset ... must handle more explicitly
uData.time = (0:size(outputSignal,2)-1)*dt;
uData.current = outputSignal(1+floor((sz(1)*sz(2)/100)*rand(200,1)),:);
plot(uData.time, uData.current);

% title('Output current');
xlabel('Time (sec)');
ylabel('pA');

end
