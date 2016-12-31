function [uData, hf] = plot(obj, pType, varargin)
% Plot function for coneMosaic base class
%
%    [uData, hf] = coneMosaic.plot(pType, varargin)
%
% There is a specialized plot() for the coneMosaicHex class that calls this
% function.
%
% Inputs:
%   pType - string, type of plot
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
%   'cone spectral qe'     - Cone pigment and macular
%
%  Pigments
%   'macular transmittance'- Graph
%   'macular absorptance'  - Graph
%   'eye spectral qe'      - Cone pigment with macular and lens
%
%  Absorptions
%   'mean absorptions'     - Image of the mean absorptions
%   'absorptions'          - Movie of the absorptions on cone mosaic
%   'movie absorptions'    - Gray scale movie of absorptions
%
%  Eye movements
%   'eye movement path'    - eye movement
%
%  Current
%   'current'              - Current movie on cone mosaic
%   'mean current'         - Image of the mean current
%   'current timeseries'   - Cone photocurrent graphs
%   'impulse response'     - Current impulse response for this integration time
%   'movie current'        - Gray scale movie of current
%
% When you present 'os ' or 'outersegment ' then we pass the arguments
% along to os.plot()
%
% Example:
%    rgc.mosaic{1}.plot(type)
%
% HJ/BW, ISETBIO TEAM, 2016

%% parse input
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('obj');
p.addRequired('pType', @isstr);               % Type of plot

p.addParameter('hf', []);                    % figure handle
p.addParameter('oi',[],@isstruct);           % Used for spectral qe

p.parse(obj,pType, varargin{:});
hf = p.Results.hf;
oi = p.Results.oi;   % Used in plotGraphs routine

uData = [];

% Find a cleaner way to check and send to os.plot().  Maybe create a parse
% argument string as in ISET.
if (length(pType) > 3 && strcmp(pType(1:3),'os '))
    obj.os.plot(pType(4:end),varargin{:});
    return;
elseif (length(pType) > 13 && strcmp(pType(1:13),'outersegment '))
    obj.os.plot(pType(14:end),varargin{:});
    return;
end

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

%% Could simplify this big switch

switch ieParamFormat(pType)
    
    % ----   Images
    case 'conemosaic'
        [uData.support, uData.spread, uData.delta, uData.mosaicImage] = ...
            conePlot(obj.coneLocs * 1e6, obj.pattern);
        imagesc(uData.mosaicImage); axis off; axis image;
    case 'meanabsorptions'
        % title('Mean number of absorptions');
        if isempty(obj.absorptions), error('no absorption data'); end
        
        % Show the data, with the gamma from the window.
        uData = mean(obj.absorptions,3);
        gdata = guidata(obj.hdl);
        gam = str2double(get(gdata.editGam,'string'));
        imagesc(uData.^gam); axis off;
        
        % Preserve the tick labels in real photons
        colormap(gray);  % Shows a numerical value
        cbar = colorbar;
        photons = str2double(get(cbar,'TickLabels')).^(1/gam);
        photons = num2str(round(photons)); set(cbar,'TickLabels',photons);
        axis image;
        
    case 'movieabsorptions'
        % Movie in gray scale
        
        % Could become cla return
        if isempty(obj.absorptions), error('no absorption data'); end
        
        % Additional arguments may be the video file name, step, and
        % FrameRate
        uData = ieMovie(obj.absorptions,varargin{:});
    case 'meancurrent'
        if isempty(obj.current), error('no photocurrent data'); end
        uData = mean(obj.current, 3);
        
        % Apply gamma.  The current is always negative.
        gdata = guidata(obj.hdl);
        gam = str2double(get(gdata.editGam,'string'));
        if max(uData(:)) > 0
            warning('Gamma correction in display is not correct');
        end
        
        % Carry on assuming current is negative pA.
        % uData = -1*(abs(uData).^gam);
        uData = abs(uData);
        if ~isequal(hf, 'none'), imagesc(uData.^gam); end
        
        axis off; colormap(flipud(gray));  % Shows a numerical value
        cbar = colorbar;
        current = -1*(abs(str2double(get(cbar,'TickLabels')).^(1/gam)));
        current = num2str(round(current)); set(cbar,'TickLabels',current);
        axis image;
    case 'moviecurrent'
        % Current movie in gray scale
        if isempty(obj.current)
            if isempty(p.Results.hf), close(hf); end
            error('no current data');
        end
        
        % Additional arguments may be the video file name, step, and
        % FrameRate
        disp('No gamma applied');
        uData = ieMovie(obj.current,varargin{:});
        
     % RGB movies - Could be resurrected some day with a different name
     %
     %     case 'absorptions'
     %         % Movie of the absorptions on the cone mosaic
     %         if isempty(obj.absorptions)
     %             % Could be come cla; return
     %             error('no absorption data');
     %         end
     %         uData = coneImageActivity(obj, hf, varargin{:});
     %     case {'current', 'photocurrent'}
     %         % Photo current movie on colored cone mosaic
     %         if isempty(obj.current)
     %             if isempty(p.Results.hf), close(hf); end
     %             error('no photocurrent data');
     %         end
     %         uData = coneImageActivity(obj, hf, 'dataType', ...
     %             'photocurrent', varargin{:});
     
        % ------ Graphs
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


% Would be used in commented out cases above, some day.
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
