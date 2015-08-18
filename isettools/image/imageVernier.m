function [Img,params] = imageVernier(params)
% Create an RGB image of a vernier line-pair
%
%   [Img, params] = imageVernier(params)
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

%% Initialize parameters
if notDefined('params'), params = []; end
if isfield(params, 'sceneSz'), sz = params.sceneSz; else sz = 64; end
if isfield(params, 'barWidth'), width = params.barWidth; else width = 1; end
if isfield(params, 'offset'), offset = params.offset; else offset = 1; end

if isfield(params, 'barColor')
    barColor = params.barColor;
else
    barColor = 0.99;
end
if isscalar(barColor), barColor = repmat(barColor, [1 3]); end

if isfield(params, 'bgColor'), bgColor=params.bgColor; else bgColor=0; end
assert(isscalar(bgColor), 'bgColor should be a scalar');

%% Create image to be shown on display
%  create 1d pattern
if isfield(params, 'pattern')
    pattern = params.pattern;
else
    if isscalar(sz), sz = [sz sz]; end
    pattern = bgColor * ones(1, sz(2));
    pattern(round((sz(2)-width)/2):round((sz(2)+width)/2)) = 1;
end
Img = image1d(pattern, 'rgb', barColor, 'rows', sz(1));

% Shift for offset
Img(1:round(end/2),:,:) = ...
    circshift(Img(1:round(end/2),:,:), [0 offset 0]);

if nargout == 2
    params.sceneSz   = sz;
    params.barWidth  = width;
    params.offset    = offset;
    params.barColor  = barColor;
    params.bgColor   = bgColor;
    params.pattern   = pattern;
end

end