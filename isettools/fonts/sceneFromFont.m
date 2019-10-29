function scene = sceneFromFont(font, display, scene, oSample, varargin)
% Create a scene from a font and display
%
% Syntax:
%   scene = sceneFromFont(...
%       [font], [display], [scene], [oSample], [varargin])
%
% Description:
%    Create a scene for a font on a display. None of the inputs are
%    required, however they allow the user to specify details not
%    necessarily available with the defaults.
%
%  Inputs:
%    font     - (Optional) Struct. A font structure.
%               Default will call fontCreate
%    display  - (Optional) Struct. A display structure.
%               Default will call displayCreate
%    scene    - (Optional) Struct. A scene structure.
%               Default will call sceneCreate
%    oSample  - (Optional) 1x2 Vector. An up-sampling rate.
%               Default is [20 20]
%    varargin - (Optional) Additional parameters, could include:
%               varargin{1} - pad size for the font bitmap (default: [])
%               varargin{2} - pad value for the font bitmap (default: [])
%
% Outputs:
%    scene    - A scene structure containing the desired font & display.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/14  BW/HJ  Vistasoft group, 2014
%    02/28/18  jnm    Formatting

%% Input arguments
if notDefined('font'),    font = fontCreate; end
if notDefined('display'), display = displayCreate('LCD-Apple'); end
if notDefined('scene'),   scene = sceneCreate; end
if notDefined('oSample'), oSample = [20 20]; end
if ~isempty(varargin),    padsz = varargin{1}; else, padsz = []; end
if length(varargin) > 1,  padval = varargin{2}; else, padval = []; end

% Initialize the display to match the scene and font properties
if ischar(display), display = displayCreate(display); end
display = displaySet(display, 'wave', sceneGet(scene, 'wave'));
if displayGet(display, 'dpi') ~= fontGet(font, 'dpi')
    warning('Adjusting display dpi to match font');
    display = displaySet(display, 'dpi', fontGet(font, 'dpi'));
end

%% Compute the high resolution display image
paddedBitmap = fontGet(font, 'padded bitmap', padsz, padval);
% vcNewGraphWin; imagesc(paddedBitmap);

np = displayGet(display, 'n primaries');
paddedBitmap = padarray(paddedBitmap, ...
    [0 0 np - size(paddedBitmap, 3)], 'post');
dRGB = displayCompute(display, paddedBitmap, oSample);
[dRGB, r, c] = RGB2XWFormat(dRGB);
spd = displayGet(display, 'spd');
wave = displayGet(display, 'wave');

% Convert the display radiance (energy) to photons and place in scene
energy = dRGB * spd';
energy = XW2RGBFormat(energy, r, c);
p = Energy2Quanta(wave, energy);
scene = sceneSet(scene, 'photons', p);   % Compressed photons

% sceneWindow(scene);

%% Adjust the scene to match the display resolution

% Adjust mean luminance to maximum Y value of display, but corrected for
% number of black pixels
wp = displayGet(display, 'white point');
nPixels = numel(paddedBitmap(:, :, 1));
p = paddedBitmap(:, :, 2);
s = sum(p(:)) / nPixels;
scene = sceneAdjustLuminance(scene, wp(2) * s);

dist = 0.5;
scene = sceneSet(scene, 'distance', dist);

% Calculate scene width in meters.  Each bitmap is 1 pixel.
dpi = displayGet(display, 'dpi');
mPerDot = dpi2mperdot(dpi, 'meters');
nDots = size(paddedBitmap, 2);
wMeters = mPerDot * nDots;
fov = atan2d(wMeters, dist);
scene = sceneSet(scene, 'fov', fov);

% Name it
scene = sceneSet(scene, 'name', fontGet(font, 'name'));

end
