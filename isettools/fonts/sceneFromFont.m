function scene = sceneFromFont(font,display,scene)
% Create a scene from a font and display
%
%  scene = sceneFromFont(font,[display='LCD-Apple'],[scene])
%
% (BW) Vistasoft group, 2014

%% Input arguments
if notDefined('font'),    font = fontCreate; end
if notDefined('display'), display = displayCreate('LCD-Apple'); end
if notDefined('scene'),   scene = sceneCreate; end

% Initialize the display to match the scene and font properties
if ischar(display), display = displayCreate(display); end
display = displaySet(display,'wave',sceneGet(scene,'wave'));
if displayGet(display,'dpi') ~= fontGet(font,'dpi')
    warning('Adjusting display dpi to match font');
    display = displaySet(display,'dpi',fontGet(font,'dpi'));
end

%% Compute the high resolution display image
paddedBitmap = fontGet(font,'padded bitmap');
np = displayGet(display, 'n primaries');
paddedBitmap = padarray(paddedBitmap, ...
                    [0 0 np - size(paddedBitmap, 3)], 'post');
dRGB       = displayCompute(display,paddedBitmap);
[dRGB,r,c] = RGB2XWFormat(dRGB);
spd  = displayGet(display,'spd');
wave = displayGet(display,'wave');

% Convert the display radiance (energy) to photons and place in scene
energy = dRGB*spd';
energy = XW2RGBFormat(energy,r,c);
p = Energy2Quanta(wave,energy);
scene = sceneSet(scene, 'cphotons', p);   % Compressed photons

% vcAddObject(scene); sceneWindow;

%% Adjust the scene to match the display resolution

% Adjust mean luminance to maximum Y value of display, but corrected
% for number of black pixels
wp = displayGet(display,'white point');
nPixels = numel(paddedBitmap(:,:,1));
p = paddedBitmap(:,:,2); s = sum(p(:))/nPixels;
scene = sceneAdjustLuminance(scene,wp(2)*s);

dist = 0.5;
scene = sceneSet(scene,'distance',dist);

% Calculate scene width in meters.  Each bitmap is 1 pixel.
dpi     = displayGet(display,'dpi');
mPerDot = dpi2mperdot(dpi,'meters');
nDots   = size(paddedBitmap,2);
wMeters = mPerDot*nDots;
fov     = atan2d(wMeters,dist);
scene   = sceneSet(scene,'fov',fov);

% Name it
scene = sceneSet(scene,'name',fontGet(font,'name'));

end
