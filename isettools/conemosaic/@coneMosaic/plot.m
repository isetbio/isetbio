function [hf, uData] = plot(obj, type, varargin)
% Plot function for coneMosaic
%    [hf, uData] = coneMosaic.plot(type, varargin)
%
% Inputs:
%   type - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'newFig' - whether to create a new figure
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
%   'absorptions'          - Movie of the absorptions
%   'current'              - Image of the mean current .... NYI
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
p.addParameter('newFig', true, @islogical);  % Generate a new figure
p.addParameter('oi',[],@isstruct);           % Used for spectral qe

p.parse(type, varargin{:});
newFig = p.Results.newFig;
oi = p.Results.oi;

uData = [];

% plot
if newFig, hf = vcNewGraphWin; else hf = []; end

switch ieParamFormat(type)
    case 'conemosaic'
        [~,~,~,coneMosaicImage] = conePlot(obj.coneLocs*1e6, obj.pattern);
        imagesc(coneMosaicImage); axis off;
    case 'conefundamentals'
        % The cone absorptance without macular pigment or lens
        c = [1 .2 0; .2 1 0; 0 0 1];
        for ii=1:3
            hold on
            plot(obj.wave, obj.cone.absorptance(:,ii), 'LineWidth', 2,'Color',c(ii,:));
        end
        grid on;
        xlabel('Wavelength (nm)'); ylabel('Cone quanta absorptance');
    case 'maculartransmittance'
        plot(obj.wave, obj.macular.transmittance, 'LineWidth', 2);
        xlabel('Wavelength (nm)'); ylabel('Macular transmittance'); grid on;
    case 'macularabsorptance'
        plot(obj.wave, obj.macular.absorptance, 'LineWidth', 2);
        xlabel('Wavelength (nm)'); ylabel('Macular absorptance'); grid on;
    case 'conespectralqe'
        % Quantum efficiency of macular pigment and cone photopigments
        plot(obj.wave, obj.qe, 'LineWidth', 2); grid on;
        xlabel('Wavelength (nm)'); ylabel('Cone quanta efficiency');
    case 'eyespectralqe'
        % Includes lens, macular pigment, and cone photopigment properties
        if isempty(oi), error('oi required for spectral qe'); end
        lensTransmittance = oiGet(oi,'lens transmittance','wave',obj.wave);
        qe = bsxfun(@times,lensTransmittance,obj.qe);
        % There is a method of setting the color order for a plot.  But we
        % didn't get it to work so we put this hack in.  We want a color
        % order set at the top. (BW/HJ)
        c = [1 .2 0; .2 1 0; 0 0 1];
        for ii=1:3
            hold on
            plot(obj.wave, qe(:,ii), 'LineWidth', 2,'Color',c(ii,:));
        end
        grid on;
        xlabel('Wavelength (nm)'); ylabel('Eye quanta efficiency');
    case 'meanabsorptions'
        if isempty(obj.absorptions)
            error('no absorption data computed');
        end
        imagesc(mean(obj.absorptions,3)); axis off;
        colorbar; title('Mean number of absorptions');
    case 'absorptions'
        % Movie of the absorptions
        if isempty(obj.absorptions)
            error('no absorption data computed');
        end
        uData.mov = coneImageActivity(obj,newFig,varargin{:});
    case 'current'
        if isempty(obj.current)
            error('no current data computed');
        end
        imagesc(obj.current(:, :, 1)); axis off;
    otherwise
        error('unsupported plot type');
end

end

%%
function mov = coneImageActivity(cMosaic,dFlag, varargin)
% Make a movie or a single image of cone absorptions on a colored mosaic
%
%   mov = coneImageActivity(cones,varargin)
% 
% cones:  sensor object (required)
% step:   How many cone absorptions to step over.  By default, if < 100
%         cone time steps (step = 1).  In general, try to show about 100
%         images. 
% dFlag:  Either a struct or a boolean
%         If a struct, dFlag contains the movie parameters.  The movie is
%         saved based on the name field in these parameters. The parameters
%         are .vname (video file name) and .FrameRate (video frame rate).
%         If dFLag is a boolean, the value indicates whether to show the
%         movie (true) or not (false).
%
% Example:
%   cones = sensorCreate; ...
%   m = coneImageActivity(cones,'step',step,'dFlag',dFlag);
%
% See Also:  conePlot, conePlotHelper, t_rgc, t_rgcIntroduction
%
% (BW) ISETBIO Team, Copyright 2015

%% Parameters

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('cMosaic',@(x) (isa(x,'coneMosaic')));
p.addRequired('dFlag',@(x) (islogical(x) || isstruct(x)));
p.addParameter('step',[],@isnumeric);
p.addParameter('coneSize',6,@isnumeric);

% Parse results
p.parse(cMosaic,dFlag,varargin{:});
cMosaic  = p.Results.cMosaic;
step     = p.Results.step;
dFlag    = p.Results.dFlag;
coneSize = p.Results.coneSize;

absorptions = cMosaic.absorptions;
if isempty(absorptions), error('No absorptions.'); end
nframes = size(absorptions,3);

% If less than 100 frames, show them all.  Otherwise show 100 frames.
if isempty(step), step = max(1,nframes/100); end

% conesToPlot = cMosaic.mosaicSize; % max(sensorGet(cMosaic,'size'));

%% Start the making the iamges

%[xy,coneType, support,spread,delta] = conePlotHelper(cMosaic, conesToPlot, coneSize);

% separation(1) = cMosaic.cone.width*10^6;   % In microns
% separation(2) = cMosaic.cone.height*10^6; 

xy = cMosaic.coneLocs*10^6; % sensorGet(sensor, 'cone xy');
% xCoords = squeeze(xy(:,1));
% yCoords = squeeze(xy(:,2));

% These are the cones we will plot.  They are chosen from the center of the
% image.
% selectConeIndices = find( ...
%     (abs(xCoords) < separation(1)*conesToPlot(2)/2) & ...
%     (abs(yCoords) < separation(2)*conesToPlot(1)/2) ...
%     );

coneType = cMosaic.pattern;    % coneType = sensorGet(sensor, 'coneType');

spread        = 0.35*coneSize; 
support       = round(spread*4*[1 1]); 
delta         = 1.5/coneSize;
% coneSubMosaic = coneType(selectConeIndices);
% xySubMosaic   = xy(selectConeIndices,:);

%
whiteBackground = false;
[~,~,~, coneMosaicImage] = conePlot(xy,coneType,support,spread,delta, whiteBackground);

% Prepare the output movie data
theseFrames = 1:step:nframes;
sz = [size(coneMosaicImage),length(theseFrames)];
mov = zeros(sz);

% Convert the image to the size of the mosaic and blur the result
wbar = waitbar(0,'creating cone movie');
g = fspecial('gaussian',6,2);      % Blurring parameters
cnt = 1;
for ii=theseFrames
    waitbar(ii/nframes,wbar);
    d = absorptions(:,:,ii);
    
    % Expand the data to the size of the coneMosaicImage
    fgrid = ffndgrid(xy,d(:),delta);
    fgrid = full(fgrid);
    
    % Scale the cone mosaic image by the adapted data
    for jj=1:3
        mov(:,:,jj,cnt) = convolvecirc(coneMosaicImage(:,:,jj) .* fgrid,g);
    end
    cnt = cnt+1;
end
close(wbar);

% Scale and gamma correct mov
% mov = ieScale(mov,0,1);
mov = mov/max(mov(:));
mov = mov .^ 0.3;
% vcNewGraphWin; imagescRGB(mov);

%% Produce the movie

if isstruct(dFlag)
    % When dFlag is a struct, show the move and save it in a file
    % Normalize tmp for visualization
    % Record a movie named by the string in showFlag
    vObj = VideoWriter(dFlag.vname);
    vObj.FrameRate = dFlag.FrameRate;
    open(vObj);
    for ii = 1:size(mov,4)
        image(mov(:,:,:,ii))
        drawnow;
        F = getframe;
        writeVideo(vObj,F);
    end
    close(vObj);
elseif dFlag
    % If it is a boolean and true, must show the movie on the screen
    for ii=1:size(mov,4)
        imshow(mov(:,:,:,ii)); 
        drawnow;
    end
else
    % dFlag must have been false.  So, we return the data but don't show it
    % or save it.
end
    
end
