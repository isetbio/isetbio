function tmp = coneImageActivity(cones,data,step,dFlag)
% Make a movie of absorptions on a color version of the cone mosaic
%
%   tmp = coneImageActivity(data,cones,step,showFlag)
%
% Example:
%
% See Also
%
% (BW) ISETBIO Team, Copyright 2015

if notDefined('data'), data = sensorGet(cones,'photons'); end
if notDefined('step'), step = max(1,size(data,3)/100); end
if notDefined('dFlag'), dFlag = false; end

conesToPlot = max(sensorGet(cones,'size'));
coneSize = 6;

[xy,coneType, support,spread,delta] = conePlotHelper(cones, conesToPlot, coneSize);
whiteBackground = false;
[~,~,~, coneMosaicImage] = conePlot(xy,coneType, support,spread,delta, whiteBackground);
% vcNewGraphWin
% imshow(coneMosaicImage); truesize;
nframes = size(data,3);

% Prepare the output movie data
theseFrames = 1:step:nframes;
sz = [size(coneMosaicImage),length(theseFrames)];
tmp = zeros(sz);

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
        tmp(:,:,jj,cnt) = convolvecirc(coneMosaicImage(:,:,jj) .* fgrid,g);
    end
    cnt = cnt+1;
end
close(wbar);

% Scale and gamma correct tmp
tmp = tmp/max(tmp(:));
tmp = tmp .^ 0.3;
% vcNewGraphWin; imagescRGB(tmp);

% Show the movie and return it (dFlag == true) Or save it as a movie file
% in based on the parameters in the dFlag struct
if isstruct(dFlag)
    % Normalize tmp for visualization
    vcNewGraphWin;
    nframes = size(tmp,4);
    % Record a movie named by the string in showFlag
    vObj = VideoWriter(dFlag.vname);
    vObj.FrameRate = dFlag.FrameRate;
    open(vObj);
    for ii = 1:step:nframes
        image(tmp(:,:,:,ii))
        drawnow;
        F = getframe;
        writeVideo(vObj,F);
    end
    close(vObj);
elseif dFlag
    % Just show the movie on the screen
    for ii=1:step:nframes
        imshow(tmp(:,:,:,ii)); 
        drawnow;
    end
end
    
end