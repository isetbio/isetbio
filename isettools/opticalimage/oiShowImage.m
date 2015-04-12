function oiShowImage(oi,displayFlag,gam)
%Render an image of the scene data
%
%    oiShowImage(oi,displayFlag,gam)
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
photons = oi.data.photons;  % Not using an oiGet() here saves memory
sz      = oiGet(oi,'size');

if isempty(photons)
    cla
    sprintf('ISET Warning:  Data are not available');
    return;
end
    
% The displayFlag flag determines how imageSPD converts the data into a
% displayed image. 
imageSPD(photons,wList,gam,sz(1),sz(2),displayFlag);
axis image; axis off

end