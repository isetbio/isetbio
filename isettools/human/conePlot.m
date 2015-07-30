function [support, spread, delta, coneMosaicImage] = ...
    conePlot(xy,coneType, support, spread, delta, whiteBackground,absorptions)
% Plot a saturated image of the human cone mosaic
%
%    [support, spread, delta, [coneMosaicImage]] = conePlot(xy,coneType, [support], [spread], [delta], [whiteBackground],[absorptions])
%
% INPUTS
%  xy:        cone xy positions (microns)
%  coneTYpe:  an integer from 1:4 where 1 means no cone (K), 2:4 are L,M,S
%            support, spread, delta are gaussian blurring parameters for
%            creating the image 
%  support:   spatial support for each cone (computed by default)
%  spread:    spatial spread for each cone  (computed by default)
%  delta:     spacing between the cones (microns)
%  whiteBackground: boolean indicating whether to generate an image with
%                  white or black background
% OUTPUTS
%  support, spread, delta - the values used to create the image
%  coneMosaicImage        - the image  
%     If this fourth output argument is present, the function returns the
%     RGB image and does NOT plotting the image
%
% These images can be compared with the ones measured by Hofer et al. in
% the J. Neuroscience paper and then published by Williams in a JOV paper.  
% I downloaded those from
%
%   http://www.journalofvision.org/5/5/5/article.aspx#tab1
%
% Those data and the plotting routine for them are in the repository under
% cone/data/williams.
%
% See also:  sensorConePlot, humanConeMosaic
%
% Example:
%  [sensor, xy, coneType] = sensorCreateConeMosaic;
%  conePlot(xy,coneType);
%
% Copyright ImagEval LLLC, 2009


if notDefined('delta'), delta = 0.4; end  % Sampling in microns
% support and spread are adjusted below, after the grid is built

% Grid the xy cone positions to the delta spacing using a fast method
fgrid = ffndgrid(xy,coneType,delta);
fgrid = full(fgrid);

% Grid the cone absorption rates the same way
if exist('absorptions','var')
    % This is the interpolation of the absorption data
    dgrid = ffndgrid(xy,absorptions(:),delta);
    dgrid = fullgrid(dgrid);
end
% Could have an else dgrid = ones(size(xy,1),1) here and then eliminate the
% else below.

% Find the positions of the empty (K) and three cone types
K = find(fgrid == 1); L = find(fgrid == 2); 
M = find(fgrid == 3); S = find(fgrid == 4);

% Wherever there is a red cone, we put (1,0,0) and so forth for the
% other cone types.
% We start with a coneImage that is nCones x 3
coneImage = zeros(numel(fgrid),3);

if exist('dgrid','var')
    coneImage(K,:) = repmat([0,0,0],length(K),1);
    coneImage(L,:) = repmat([dgrid(L),0,0],length(L),1);
    coneImage(M,:) = repmat([0,dgrid(M),0],length(M),1);
    coneImage(S,:) = repmat([0,0,dgrid(S)],length(S),1);
else
    % We set the L cone rows to red, ... and so forth
    coneImage(K,:) = repmat([0,0,0],length(K),1);
    coneImage(L,:) = repmat([1,0,0],length(L),1);
    coneImage(M,:) = repmat([0,1,0],length(M),1);
    coneImage(S,:) = repmat([0,0,1],length(S),1);
end

% Reshape the image to its actual (row,col,3) size
coneImage = reshape(coneImage,size(fgrid,1),size(fgrid,2),3);
% mp = [0 0 0 ; 1 0 0 ; 0 1 0; 0 0 1]; image(fgrid); colormap(mp)

% Blur the image by a Gaussian - we set blur and support here.
if notDefined('spread')
    l = find(fgrid(1,:) > 0);  % Digital image spacing
    spread = (l(2)-l(1))/3;
end
if notDefined('support'), support = round(3*[spread spread]); end


if notDefined('whiteBackground'), whiteBackground = false; end
if (whiteBackground)
    g = fspecial('gaussian',support,spread*0.87);
    g = g/max(g(:));
    g(g<0.1) = 0;
    g = 1.5*g.^0.3;
    g(g>1) = 1.0;
else
    g = fspecial('gaussian',support,spread);
end

tmp = imfilter(coneImage,g);

if (whiteBackground)
    tmp = tmp/max(tmp(:));
    indices = all(tmp < .75,3);
    tmp(repmat(indices,[1,1,3])) = 1;
end

if (nargout < 4)
    % Show the image
    h = vcNewGraphWin;
    set(h,'Name','ISET: Human cone mosaic');
    imagescRGB(tmp);
else
    % return the image
    coneMosaicImage = tmp/max(tmp(:));
end