function mov = coneImageActivity(cones,varargin)
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
p.addRequired('cones');
p.addParameter('step',[],@isnumeric);
p.addParameter('dFlag',false,@(x) (islogical(x) || isstruct(x)));
p.addParameter('coneSize',6,@isnumeric);
p.parse(cones,varargin{:});

% Parse results
cones = p.Results.cones;
step  = p.Results.step;
dFlag = p.Results.dFlag;
coneSize = p.Results.coneSize;

data = sensorGet(cones,'photons'); 
if isempty(data), error('No photon absorptions in the sensor object.'); end

% If less than 100 frames, show them all.  Otherwise show 100 frames.
if isempty(step), step = max(1,size(data,3)/100); end

conesToPlot = max(sensorGet(cones,'size'));

%% Start the making the iamges
[xy,coneType, support,spread,delta] = conePlotHelper(cones, conesToPlot, coneSize);
whiteBackground = false;
[~,~,~, coneMosaicImage] = conePlot(xy,coneType, support,spread,delta, whiteBackground);

% vcNewGraphWin
% imshow(coneMosaicImage); truesize;
nframes = size(data,3);

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
    d = data(:,:,ii);
    
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
    vcNewGraphWin;
    nframes = size(mov,4);
    % Record a movie named by the string in showFlag
    vObj = VideoWriter(dFlag.vname);
    vObj.FrameRate = dFlag.FrameRate;
    open(vObj);
    for ii = 1:step:nframes
        image(mov(:,:,:,ii))
        drawnow;
        F = getframe;
        writeVideo(vObj,F);
    end
    close(vObj);
elseif dFlag
    % If it is a boolean and true, uust show the movie on the screen
    for ii=1:step:nframes
        imshow(mov(:,:,:,ii)); 
        drawnow;
    end
else
    % dFlag must have been false.  So, we return the data but don't show it
    % or save it.
end
    
end