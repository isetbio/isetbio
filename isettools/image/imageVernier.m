function [I, params] = imageVernier(params, varargin)
% Create an RGB image of a vernier line-pair
%
%   [I, params] = imageVernier(params)
%
% The image, typically a pair of lines that are offset is created from the
% parameters.  It is possible, however, to simply send in a pattern.
%
% The default value for vernier image params are
%
%   p.sceneSz   = 64    - Scene size in pixels
%   p.barWidth  = 1     - Bar width in pixels
%   p.barLength = p.sceneSz(1)
%   p.barColor  = 1     - Bar color (either scalar or RGB)
%   p.gap       = 0     - Pixel gap between the upper and lower images
%
%   p.offset    = 1     - Offset between the top and bottom pattern
%   p.bgColor   = 0     - Background color (either scalar or RGB)
%   p.pattern   = []    - A 1D pattern
%
% If the pattern is not specified, the program creates vernier image based
% on the basic one line pattern. If pattern is specified, the image is
% created based on pattern and some parameters might not be used
%   
% Examples:
% You can run these in series to see the effect of each parameter
%
%   p = vernierP;    img = imageVernier(p); imshow(img)
%   p.sceneSz = 256; img = imageVernier(p); imshow(img);
%   p.gap  = 6;      img = imageVernier(p); imshow(img);
%   p.bgColor = 0.2; img = imageVernier(p); imshow(img);
%   p.barLength = 16; img = imageVernier(p); imshow(img);
%   p.barColor = [1 0 0]; img = imageVernier(p); imshow(img);
%   p.bgColor = [0 1 1]; img = imageVernier(p); imshow(img);
%
%   s = sceneCreate('vernier','display',p);
%   ieAddObject(s); sceneWindow;
%
%   p.display = displayCreate('OLED-Sony','dpi',300);
%   s = sceneCreate('vernier','display',p);
%   ieAddObject(s); sceneWindow;
%
%   x = (-63:64)/128; f = 2;
%   p.pattern = 0.5*cos(2*pi*f*x) + 0.5;
%   img = imageVernier(p); imshow(img);
%
% HJ/BW ISETBIO Team Copyright 2015

%% Parse input parameters
p = inputParser; p.KeepUnmatched = true;

p.addParameter('sceneSz', 64, @(x) isnumeric(x));
p.addParameter('barWidth', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('barLength', [], @isnumeric);

p.addParameter('offset', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('gap', 0, @isscalar);   % Pixel gap between upper and lower

p.addParameter('barColor', 1, @isnumeric);
p.addParameter('bgColor', 0, @(x) isnumeric(x));

% A 1D pattern to use instead of the default line
p.addParameter('pattern', []);

p.parse(params, varargin{:});
params = p.Results;
sz     = params.sceneSz;
width  = params.barWidth;
offset = params.offset;
barColor = params.barColor;
bgColor  = params.bgColor;
barLen   = params.barLength;
gap      = params.gap;

if isempty(barLen), params.barLength = sz(1); barLen = sz(1); end
if isscalar(barColor), barColor = repmat(barColor, [1 3]); end
if isscalar(bgColor),  bgColor = repmat(bgColor, [1 3]); end

%% Create 1d pattern
if ~isempty(params.pattern)
    pattern = params.pattern;
    if ismatrix(pattern), pattern = repmat(pattern, [1 1 3]); end
else
    % We are building a 1D pattern that is a line with the appropriate
    % width, barColor and background color.
    if isscalar(sz), sz = [sz sz]; end
    pattern = bsxfun(@times, reshape(bgColor,1,1,[]), ones(1,sz(2),3));
    barIndx = round((sz(2)-width)/2):round((sz(2)-width)/2+width-1);
    for ii = 1 : 3, pattern(1, barIndx, ii) = barColor(ii); end
end

% Convert the 1D pattern into an image of the appropriate length
I = repmat(pattern, [barLen 1 1]);

% Insert the gap dealing with even and odd considerations
I = insertGap(I,gap);

% Shift the upper rows by the offset
I(1:round(end/2),:,:) = circshift(I(1:round(end/2),:,:), [0 offset 0]);

% Pad rows with background color to get desired image size
I = padarray(I, [ceil((sz(1)-barLen)/2) 0 0], nan);
I = I(1:sz(1), :, :);
for ii = 1 : 3
    curImg = I(:, :, ii);
    curImg(isnan(curImg)) = bgColor(ii);
    I(:, :, ii) = curImg;
end

end

%---------
function I = insertGap(I,gap)
% Insert a gap into the image by placing NaNs into the rows
%
% gap is the number of rows
% The even and odd cases are handled separately

sz = size(I);
if isodd(sz(1)) && isodd(gap)
    % We put NaNs in the middle row and the appropriate number above and
    % below that row 
    mid = (sz(1) + 1)/2;
    rows = ((1:gap) - (gap+1)/2) + mid;    
elseif ~isodd(sz(1)) && ~isodd(gap)
    mid = sz(1)/2;
    rows = (1:gap) - gap/2 + mid; 
else
    error('Bad row size %d, gap size %d pair.',sz(1),gap)
end

% Put Nans in the gap rows.  They will get filled in with the bgColor 
for ii=1:3
    I(rows,:,ii) = nan;
end

end