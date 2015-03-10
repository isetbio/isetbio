function p = plotChromaticityDiagram(background,imageLocation)
% Doesn't seem to work - use chromaticityPlot.
%
%    Draws the xy chromaticity diagram on a new figure
%
%    p = plotChromaticityDiagram(background,imageLocation)
%
% INPUTS"
%  background - is the background color; can only be be 'white' or 'black'. 
%  These are the most convenient 2 options.  Matlab figure properties 
%  (grids, colors, etc.) can be changed after the chromaticity diagram
%  is drawn.
%
%  imageLocation - is the location (fullfile) of a previously saved image.
%  If this is not specified, we will look in the folder where this function
%  resides. If we don't find the image, we will create and save an image
%  for the given background color. In subsequent runs, the routine will not
%  create a new diagram; it will simply load the previously saved image
%  (and work much faster).
%
% Examples:
%  background = 'white'
%  fig_handle = plotChromaticityDiagram(background) 
%  fig_handle = plotChromaticityDiagram(background,imageLocation)
%
% Copyright ImagEval 2011

error('See chromaticityPlot')

return

%%

if notDefined('imageLocation')
    imageLocation = fileparts(mfilename('fullpath'));
end

if notDefined('background')
    background = 'white';
elseif ~strcmp(background,'black') && ~strcmp(background,'white')
    disp('I am going to draw the diagram on a white background');
    background = 'white';
end

filename = fullfile(imageLocation,sprintf('chromaticityDiagram_%s.mat',background));

% If such a file exists, we will use it to draw the diagram
if exist(filename, 'file') == 2
    p = plotChromaticityDiagram(filename, background);
    return
end
% Else, we create a new diagram and save it at imageLocation

% Set resolution of image. Bigger is better but slower.
n_pixels  = 512;
x = linspace(0.001,1,n_pixels); % avoid divide by 0 by skipping the origin
y = linspace(0.001,1,n_pixels);

% Set a value for Y
Y_val = 40;

% get xy coordinates and the appropriate set of xyY points
[xx,yy] = meshgrid(x,y);
[n_rows,n_cols] = size(xx);

n_points = n_rows*n_cols;

xy = zeros(n_rows,n_cols,2);
xy(:,:,1) = xx; xy(:,:,2) = yy;
xy = reshape(xy,n_points,2);


%%

% convert xyY to XYZ
xyY = zeros(n_points,3);
xyY(:,1:2) = xy;
xyY(:,3)   = repmat(Y_val,n_points,1);

Y = xyY(:,3);
X = (xyY(:,1) ./ xyY(:,2)) .* Y;
Z = Y .* ( ( ones(length(xyY),1) - xyY(:,1) - xyY(:,2) ) ./ xyY(:,2) );
XYZ = [X Y Z];

% convert XYZ to sRGB
matrix = [ 2.3647  -0.8965 -0.4681;...
          -0.51512  1.4264  0.0888;...
           0.0052  -0.0144  1.0092];

RGB   = (matrix*XYZ')';
m_RGB = max(RGB,[],2);
RGB   = RGB./repmat(m_RGB,1,3);
sRGB  = 255*( ieClip(RGB,0, 1) ).^2.4;

%%
% Color the points inside XYZ locus; points outside the locus are the color
% of the background

wave = 380:700;
XYZ  = ieReadSpectra('XYZ',wave);
spectrumLocus = chromaticity(XYZ);

color_me = inpolygon(xy(:,1),xy(:,2),...
    spectrumLocus(:,1),spectrumLocus(:,2));

srgb_img = 255*ones(size(sRGB));
srgb_img(color_me,:) = sRGB(color_me,:);

if strcmp(background,'black')
    srgb_img(~color_me,:) = 0;
end

srgb_img = XW2RGBFormat(srgb_img,n_rows,n_cols);

save(filename,'srgb_img','x','y');

p = plotChromaticityDiagram(filename, background);

return


function p = plotChromaticityDiagram(filename,background)
% Render the horseshoe and draw spectrum locus

% load and normalize image of xy diagram 
tmp = load(filename);
max_val = max(max(tmp.srgb_img(:,:,3)));
rgb_img = tmp.srgb_img/max_val;
x = tmp.x;
y = tmp.y;

% Set figure colors according to background
switch background
    
    case 'white'
        
        fg = [0.3 0.3 0.3];
        bg = [1 1 1];
        
    case 'black'
        
        fg = [0.7 0.7 0.7];
        bg = [0 0 0];
        
end

% Render the image

figure, 
p = imagesc(x,y,rgb_img.^0.8);
axis image;
set(gcf,'Color', bg)
set(gca,'Color', bg,'xcolor', fg, 'ycolor', fg,...
    'xgrid','on','ygrid','on','xlim',[0, 0.75],'ylim',[0,0.85],...
    'xcolor',fg, 'ycolor', fg,...
    'XMinorGrid', 'off', 'YMinorGrid', 'off',...
    'MinorGridLineStyle', ':','GridLineStyle',':',...
    'Fontsize',14,'Fontweight','normal',...
    'ydir', 'normal');
hold on

% Draw the spectrum locus

wave = 380:700;
XYZ  = ieReadSpectra('XYZ',wave);
spectrumLocus = chromaticity(XYZ);

plot(spectrumLocus(:,1),spectrumLocus(:,2),...
    'Linewidth', 1','Color', fg);
hold on
line([spectrumLocus(1,1),spectrumLocus(end,1)],...
    [spectrumLocus(1,2),spectrumLocus(end,2)],...
    'Linewidth', 1, 'Color', fg );

xlabel('x', 'Fontsize', 14, 'color', fg);
ylabel('y', 'Fontsize', 14, 'color', fg);

return