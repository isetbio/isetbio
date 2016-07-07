function img = oiShowImage(oi,displayFlag,gam)
% Render an image of the scene data
%
%    img = oiShowImage(oi,displayFlag,gam)
%
% oi: optical image
% displayFlag: <= 0 means no display. (default = 1, visible RGB and show)
%     Absolute value used to guide the type of rendering.  Thus,
%     if abs(displayFlag) = 
%              1:  Typical visible band to RGB
%              2:  Gray scale image, used for SWIR, NIR
% gam: Exponent for gamma correction
%
% Returns the display image
%
% Examples:
%   vcNewGraphWin;
%   displayFlag = 1; oiShowImage(oi,displayFlag); % RGB
%   displayFlag = 2; oiShowImage(oi,displayFlag); % Gray scale
%   displayFlag = 1; oiShowImage(oi,displayFlag, 0.5); % RGB, gamma
%
% Copyright ImagEval Consultants, LLC, 2003.

% TODO:  Shouldn't we select the axes for rendering here?  There is only
% one axis in the scene and oi window. But if we ever go to more, this
% routine should  say which axis should be used for rendering.

if isempty(oi), cla; return;  end

if notDefined('gam'), gam = 1; end
if notDefined('displayFlag'), displayFlag = 1; end

% Force to lower case and no spaces
wList   = oiGet(oi,'wavelength');
photons = [];
if checkfields(oi,'data','photons')
    photons = oi.data.photons;  % Not using an oiGet() here saves memory
end
sz      = oiGet(oi,'size');

if isempty(photons)
    handles = ieSessionGet('opticalimagehandle');
    % axes(get(handles);
    cla
    ieInWindowMessage('No spectral irradiance data available',handles,[]);
    return;
end
    
% The displayFlag flag determines how imageSPD converts the data into a
% displayed image. 
img = imageSPD(photons,wList,gam,sz(1),sz(2),displayFlag);
axis image; axis off

end