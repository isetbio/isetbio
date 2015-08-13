function [Img,params] = imageVernier(params)
% Create an RGB image of a vernier line-pair
%
%   [Img, params] = imageVernier(params)
%
% The vernier image params are
%   p.sceneSz   = 64
%   p.barWidth  = 1
%   p.offset    = 1;
%   p.lineSpace = inf
%   p.barColor  = 0.99
%   p.bgColor   = 0;
%   p.pattern   = ones(...) % For the future
%   
% Examples:
%   img = imageVernier(p); imshow(img)
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
% HJ/BW ISETBIO Team Copyright 2015
% init parameters from params


% TODO
%  We would like to put in a 1D function that could be the pattern.  RIght
%  now we only have a line.  But if we put in a little Gaussian 1D or
%  anything we could shift that instead of the line

%% Initialize parameters
if isfield(params, 'sceneSz'),   sz = params.sceneSz; else sz = 64; end
if isfield(params, 'barWidth'),  width = params.barWidth; else width = 1; end
if isfield(params, 'offset'),    offset = params.offset; else offset = 1; end

if isfield(params, 'barColor')
    barColor = params.barColor;
else
    barColor = 0.99;
end
if isscalar(barColor), barColor = repmat(barColor, [1 3]); end
if isfield(params, 'bgColor')
    bgColor = params.bgColor;
else
    bgColor = 0;
end
if isscalar(bgColor), bgColor = repmat(bgColor, [1 3]); end

%% Create image to be shown on display
if isscalar(sz), sz = [sz sz]; end
Img = repmat(reshape(bgColor,[1 1 3]),[sz 1]);

% Center column, accounting for the line spacing
cc = round(sz(2)/2);

% We change the width
width = width - 1;
for jj = 1 : length(cc)
    barCols = max(round(cc(jj)-width/2),1) : ...
        min(round(cc(jj)+width/2),sz(2));
    for ii = 1 : 3
        Img(:, barCols, ii) = barColor(ii);
    end
end

% Shift for offset
Img(1:round(end/2),:,:) = ...
    circshift(Img(1:round(end/2),:,:), [0 offset 0]);

if nargout == 2
    params.sceneSz   = sz;
    params.barWidth  = width + 1;
    params.offset    = offset;
    params.barColor  = barColor;
    params.bgColor   = bgColor;
    params.pattern   = [];   % For the future 1D pattern
end


end


