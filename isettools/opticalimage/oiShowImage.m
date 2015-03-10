function rgb = oiShowImage(oi,displayFlag,gam)
%Render an image of the scene data
%
%    rgb = oiShowImage(oi,displayFlag,gam)
%
% The rendering can be either of photons or energy values. This is called
% from the sceneWindow, so that axes are in that window.  If you call this
% from the command line, a new figure is displayed.
%
% Examples:
%   oiShowImage(oi,'photons')
%   oiShowImage(oi,'energy')
%   oiShowImage(oi,'luminance')
%
% Copyright ImagEval Consultants, LLC, 2003.

% TODO:  Shouldn't we select the axes for rendering here?  There is only
% one axis in the scene and oi window. But if we ever go to more, this
% routine should  say which axis should be used for rendering.

if isempty(oi), cla; return;  end

if notDefined('gam'), gam = 1; end
if notDefined('displayFlag'), displayFlag = 1; end

% Force to lower case and no spaces
wList = oiGet(oi,'wavelength');
img   = oiGet(oi,'photons');
sz    = oiGet(oi,'size');

if isempty(img)
    cla
    sprintf('ISET Warning:  Data are not available');
    return;
end
    
% This displays the image in the GUI.  The displayFlag flag determines how
% imageSPD converts the data into a displayed image. The data in img are
% in RGB format.
rgb = imageSPD(img,wList,gam,sz(1),sz(2),displayFlag);
axis image; axis off

end