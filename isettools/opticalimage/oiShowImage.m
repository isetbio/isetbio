function img = oiShowImage(oi, displayFlag, gam)
% Render an image of the scene data
%
% Syntax:
%   img = oiShowImage(oi, displayFlag, gam)
%
% Description:
%    Render & return the display image of the scene data.
%
%    There are examples contained in the code below. To access, type 'edit
%    oiShowImage.m' into the Command Window.
%
% Inputs:
%    oi          - (Optional) Struct. An optical image structure. If not
%                  provided, function will clear everything and return.
%    displayFlag - (Optional) Numeric. The value will indicate whether or
%                  not to display the scene, and if to display, how. The
%                  values are read using the absolute values. Default is 1.
%                  The values are as follows:
%          0: If the value <= 0, No display.
%          1: Default. Visible RGB display
%          2: Gray scale image, used for SWIR, NIR
%    gam         - (Optional). Scalar Numeroic. Exponent for gamma
%                  correction. Default is 1.
%
% Outputs:
%    img         - The created scene image
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Shouldn't we select the axes for rendering here? There is only
%      one axis in the scene and oi window. But if we ever go to more, this
%      routine should say which axis should be used for rendering.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/06/18  jnm  Formatting
%    06/24/19  JNM  Minor formatting adjustments

% Examples:
%{
    oi = oiCreate;
    scene = sceneCreate;
    oi = oiSet(oi, 'photons', oiCalculateIrradiance(scene, oi))
    vcNewGraphWin;
    displayFlag = [1 2 1];
    oiShowImage(oi, displayFlag(1)); % RGB
    vcNewGraphWin;
    oiShowImage(oi, displayFlag(2)); % Gray scale
    vcNewGraphWin;
    oiShowImage(oi, displayFlag(3), 0.5); % RGB, gamma
%}

if isempty(oi), cla; return;  end
if notDefined('gam'), gam = 1; end
if notDefined('displayFlag'), displayFlag = 1; end

% Force to lower case and no spaces
wList = oiGet(oi, 'wavelength');
photons = [];

% if we don't use an oiGet here we save memory.
if checkfields(oi, 'data', 'photons'), photons = oi.data.photons; end
sz = oiGet(oi, 'size');

if isempty(photons)
    handles = ieSessionGet('opticalimagehandle');
    % axes(get(handles);
    cla
    ieInWindowMessage(...
        'No spectral irradiance data available', handles, []);
    return;
end

% The displayFlag flag determines how imageSPD converts the data into a
% displayed image.
img = imageSPD(photons, wList, gam, sz(1), sz(2), displayFlag);
axis image;
axis off

end