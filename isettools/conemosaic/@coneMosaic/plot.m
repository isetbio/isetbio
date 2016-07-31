function [uData, hf] = plot(obj, type, varargin)
% Plot function for coneMosaic
%    [uData, hf] = coneMosaic.plot(type, varargin)
%
% Inputs:
%   type - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'hf' - figure handle or control structure, the meaning of value is
%             []            - create plot in new figure
%             figure handle - plot in figure specified by hf
%             'none'        - don't plot
%
% Outputs:
%   hf    - figure handle
%   uData - Computed data
%
% Plot type can be chosen from
%   'cone mosaic'          - Color image of the cone arrangement
%   'cone fundamentals'    - Cone pigment without macular or lens
%   'macular transmittance'- Graph
%   'macular absorptance'  - Graph
%   'cone spectral qe'     - Cone pigment and macular
%   'eye spectral qe'      - Cone pigment with macular and lens
%   'mean absorptions'     - Image of the mean absorptions
%   'absorptions'          - Movie of the absorptions on cone mosaic
%   'movie absorptions'    - Gray scale movie of absorptions
%   'mean current'         - Image of the mean current
%   'current'              - Current movie on cone mosaic
%   'movie current'        - Gray scale movie of current
%   'eye movement path'    - eye movement
%   'current timeseries'   - Cone photocurrent graphs
%
% Example:
%   cm = coneMosaic;
%   cm.plot('macular transmittance')
%   cm.plot('spectral qe','oi',oi)
%   cm.plot('absorptions','dFlag',true);  % Display the movie
%
% HJ/BW, ISETBIO TEAM, 2016


% parse input
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('type', @isstr);               % Type of plot
p.addParameter('hf', []);                    % figure handle
p.addParameter('oi',[],@isstruct);           % Used for spectral qe

p.parse(type, varargin{:});
hf = p.Results.hf;
oi = p.Results.oi;

uData = [];

% plot
if isempty(hf), hf = vcNewGraphWin;  
elseif isgraphics(hf, 'figure'), figure(hf); 
elseif isgraphics(hf, 'axes'), axes(hf);
end

% set color order so that LMS plots as RGB
if ~isequal(hf, 'none')
    co = get(gca, 'ColorOrder');
    if isgraphics(hf,'axes')
        set(get(hf,'parent'),'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :)) 
    else  % Figure
        set(hf, 'DefaultAxesColorOrder', co([2 5 1 3 4 6 7], :));
    end
end

switch ieParamFormat(type)
    case 'conemosaic'
        % we first plot to a hidden figure
        hiddenF = figure('visible', 'off');
        
        % get position of the cones in um
        x = obj.coneLocs(:, 1) * 1e6; y = obj.coneLocs(:, 2) * 1e6;
        coneWidth = obj.pigment.width * 1e6;
        coneHeight = obj.pigment.height * 1e6;
        pdWidth = obj.pigment.pdWidth * 1e6;
        
        % set image to high resolution
        curUnits = get(gca, 'Units'); set(gca, 'Units', 'Points');
        set(hiddenF, 'Position', [0 0 1200 1200]);
        pos = [100 100 1000 1000];
        
        % adjust aspect ratio
        if pos(4) < pos(3) * obj.height / obj.width;
            pos(3) = pos(4) * obj.width / obj.height;
        else
            pos(4) = pos(3) * obj.height / obj.width;
        end
        
        % compute conversion factor between meters and points
        set(gca, 'Position', pos); set(gca, 'Units', curUnits);
        meter2point = @(xm) xm * pos(3) / (max(x) - min(x) + coneWidth);
        
        % scatter plot the cones
        coneColor = [0 0 0; 0.9 0.1 0.1; 0.1 0.8 0.2; 0.1 0.1 1];
        uData.s = scatter(x, y, meter2point(pdWidth)^2, ...
            coneColor(obj.pattern, :), 'fill', 'o');
        
        % set limits of x/y axis based on the position of cones
        xlim([min(x) - coneWidth/2,  max(x) + coneWidth/2]);
        ylim([min(y) - coneHeight/2, max(y) + coneHeight/2]);
        axis equal; xlabel('Position (um)'); ylabel('Position (um)');
        set(gca, 'Color', 'black');
        set(gca, 'FontSize', 20);
        
        % get image
        set(gca, 'Units', 'Pixels'); ti = get(gca, 'TightInset');
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        set(gca, 'Units', curUnits);
        frame = getframe(gca, rect);
        uData.mosaicImage = frame.cdata;
        
        % show image if hf is given
        if isgraphics(hf, 'figure'), figure(hf); imshow(frame.cdata);
        else isgraphics(hf, 'axes'), axes(hf); imshow(frame.cdata);
        end
        
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
    case 'meanabsorptions'
        if isempty(obj.absorptions)
            if isempty(p.Results.hf), close(hf); end
            error('no absorption data');
        end
        uData = mean(obj.absorptions,3);
        if ~isequal(hf, 'none')
            imagesc(uData); axis off; colorbar; 
            title('Mean number of absorptions');
        end
        colormap(gray);  % Shows a numerical value
        axis image;
    case 'absorptions'
        % Movie of the absorptions on the cone mosaic
        if isempty(obj.absorptions)
            if isempty(p.Results.hf), close(hf); end
            error('no absorption data');
        end
        uData = coneImageActivity(obj, hf, varargin{:});
    case 'movieabsorptions'
        % Movie in gray scale
        if isempty(obj.absorptions)
            if isempty(p.Results.hf), close(hf); end
            error('no absorption data');
        end
        % Additional arguments may be the video file name, step, and
        % FrameRate
        uData = ieMovie(obj.absorptions,varargin{:});
    case 'meancurrent'
        if isempty(obj.current)
            if isempty(p.Results.hf), close(hf); end
            error('no photocurrent data computed');
        end
        uData = mean(obj.current, 3);
        if ~isequal(hf, 'none')
            imagesc(uData); axis off; colorbar;
            title('Mean photocurrent (pA)');
        end
        colormap(gray); % Shows a numerical value
        axis image;
    case {'current', 'photocurrent'}
        % Photo current movie on colored cone mosaic
        if isempty(obj.current)
            if isempty(p.Results.hf), close(hf); end
            error('no photocurrent data');
        end
        uData = coneImageActivity(obj, hf, 'dataType', ...
            'photocurrent', varargin{:});
    case 'moviecurrent'
        % Current movie in gray scale
        if isempty(obj.current)
            if isempty(p.Results.hf), close(hf); end
            error('no current data');
        end
        % Additional arguments may be the video file name, step, and
        % FrameRate
        uData = ieMovie(obj.current,varargin{:});
    case {'currenttimeseries'}
        % Photocurrent time series of selected points.
        % Need a way to choose which points!
        if isempty(obj.current)
            if isempty(p.Results.hf), close(hf); end
            error('no photocurrent data');
        end
        uData = plotCurrentTimeseries(obj,varargin{:});
        
    case {'empath', 'eyemovementpath'}
        plot(obj.emPositions(:, 1), obj.emPositions(:, 2));
        grid on; xlabel('Horizontal position (cones)');
        ylabel('Vertical position (cones)');
    otherwise
        error('unsupported plot type');
end

end

function mov = coneImageActivity(cMosaic, hf, varargin)
% Make a movie or a single image of cone absorptions on a colored mosaic
%
%   mov = coneImageActivity(coneMosaic,varargin)
% 
% cones:  coneMosaic class object=
% hf:     figure handle or 'none' or a struct
%         If a struct, hf contains the movie parameters.  The movie is
%         saved based on the name field in these parameters. The parameters
%         are .vname (video file name) and .FrameRate (video frame rate).
%         If dFLag is a boolean, the value indicates whether to show the
%         movie (true) or not (false).
%
%
% HJ/BW, ISETBIO Team, 2016

% parse inputs
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


%%
function uData = plotCurrentTimeseries(obj,varargin)
% Pull out the time series of the photo current and plot it

% Temporal samples
dt = obj.sampleTime;

% The current
outputSignalTemp = obj.current;
sz = size(outputSignalTemp);

% Reshape for plotting
outputSignal = reshape(outputSignalTemp,sz(1)*sz(2),sz(3));

% Plot a random subset ... must handle more explicitly
uData.time = (0:size(outputSignal,2)-1)*dt;
uData.current = outputSignal(1+floor((sz(1)*sz(2)/100)*rand(200,1)),:);
plot(uData.time, uData.current);

title('Output current');
xlabel('Time (sec)');
ylabel('pA');

end
