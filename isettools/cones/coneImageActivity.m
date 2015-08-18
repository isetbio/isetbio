function tmp = coneImageActivity(cones,data,step,showFlag)
% Make a move of activity on the color cone mosaic
%
%   tmp = coneImageActivity(data,cones,step,showFlag)
%
% Example:
%
% See Also
%
% (BW) ISETBIO Team, Copyright 2015

if ieNotDefined('data'), data = sensorGet(cones,'photons'); end
if ieNotDefined('step'), step = max(1,size(data,3)/100); end
if ieNotDefined('showFlag'), showFlag = false; end

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
wbar = waitbar(0,'cone movie');
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
    if showFlag, imshow(tmp(:,:,:,cnt).^0.3); drawnow; end
    cnt = cnt+1;
end
close(wbar);

end