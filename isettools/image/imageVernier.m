function [I, params] = imageVernier(params, varargin)
% Create an RGB image of a vernier line-pair
%
%   [I, params] = imageVernier(params)
%
% The default value for vernier image params are
%   p.sceneSz   = 64
%   p.barWidth  = 1
%   p.offset    = 1
%   p.lineSpace = inf
%   p.barColor  = 0.99
%   p.bgColor   = 0
%   p.pattern   = []
%
% If the pattern is not specified, the program craetes vernier image based
% on the basic one line pattern. If pattern is specified, the image is
% created based on pattern and some parameters might not be used
%   
% Examples:
%   img = imageVernier(); imshow(img)
%   p.sceneSz = 256; img = imageVernier(p); imshow(img);
%   p.bgColor = 0.5; img = imageVernier(p); imshow(img);
%
%   s = sceneCreate('vernier','display',p);
%   vcAddObject(s); sceneWindow;
%
%   p.display = displayCreate('OLED-Sony','dpi',300);
%   s = sceneCreate('vernier','display',p);
%   vcAddObject(s); sceneWindow;
%
%   x = (-63:64)/128; f = 2;
%   p.pattern = 0.5*cos(2*pi*f*x) + 0.5;
%   img = imageVernier(p); imshow(img);
%
% HJ/BW ISETBIO Team Copyright 2015

% Parse input parameters
p = inputParser; p.KeepUnmatched = true;
p.addParameter('sceneSz', 64, @(x) isnumeric(x));
p.addParameter('barWidth', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('offset', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('barLength', [], @isnumeric);
p.addParameter('barColor', 1, @isnumeric);
p.addParameter('bgColor', 0, @(x) isnumeric(x));
p.addParameter('pattern', []);

p.parse(params, varargin{:});
params = p.Results;
sz = params.sceneSz;
width = params.barWidth;
offset = params.offset;
barColor = params.barColor;
bgColor=params.bgColor;
barLen = params.barLength;

if isempty(barLen), params.barLength = sz(1); barLen = sz(1); end
if isscalar(barColor), barColor = repmat(barColor, [1 3]); end
if isscalar(bgColor), bgColor = repmat(bgColor, [1 3]); end

% Create 1d pattern
if ~isempty(params.pattern)
    pattern = params.pattern;
    if ismatrix(pattern), pattern = repmat(pattern, [1 1 3]); end
else
    if isscalar(sz), sz = [sz sz]; end
    pattern = bsxfun(@times, reshape(bgColor,1,1,[]), ones(1,sz(2),3));
    barIndx = round((sz(2)-width)/2):round((sz(2)-width)/2+width-1);
    for ii = 1 : 3, pattern(1, barIndx, ii) = barColor(ii); end
end

% Create image of the center bar part
I = repmat(pattern, [barLen 1 1]);

% Shift for offset
I(1:round(end/2),:,:) = circshift(I(1:round(end/2),:,:), [0 offset 0]);

% Pad rows with background color to get desired size
I = padarray(I, [ceil((sz(1)-barLen)/2) 0 0], nan);
I = I(1:sz(1), :, :);
for ii = 1 : 3
    curImg = I(:, :, ii);
    curImg(isnan(curImg)) = bgColor(ii);
    I(:, :, ii) = curImg;
end

end