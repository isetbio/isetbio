function outImg = render_oled_samsung(inImg, d, sz)
% render image for Samsung S-strip display
%
% Syntax:
%   I = render_oled_samsung(inImg, [d], [sz])
%
% Description:
%    This is the render function for OLED-Samsung display
%    The display can be created with
%       d = displayCreate('OLED-Samsung');
%
%    The sub-pixel design is Samsung S-strip, whose repeating unit contains
%    2x2 pixels
%
% Inputs:
%    inImg  - Matrix. An input image matrix.
%    d      - (Optional) Struct. A display structure. Default is created by
%             calling displayCreate.
%    sz     - (Optional) Numeric. The pixels per dixel. Default is to not
%             provide and instead use over sampling.
%
% Output:
%   outImg - Matrix. The rendered image.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   displayCompute, render_lcd_samsung_rgbw
%

% History:
%    XX/XX/14  HJ   ISETBIO TEAM, 2014
%    07/18/19  JNM  Formatting

% Examples:
%{
    d = displayCreate('OLED-Samsung');
    I = 0.5 * (sin(2 * pi * (1:32) / 32) + 1);
    I = repmat(I, 32, 1);
    outI = render_oled_samsung(I, d);
%}

%% Init
if notDefined('inImg'), error('input image required'); end
if notDefined('d'), d = displayCreate('OLED-Samsung'); end

if notDefined('sz')
    s = displayGet(d, 'over sample');
    controlMap = displayGet(d, 'dixel control map');
else
    s = round(sz ./ displayGet(d, 'pixels per dixel'));
    controlMap = displayGet(d, 'dixel control map', sz);
end

%% Render
% Get parameters from display structure
pixels_per_dixel = displayGet(d, 'pixels per dixel');
nprimaries = size(inImg, 3);

% process by block
outImg = zeros(size(inImg, 1) * s(1), size(inImg, 2) * s(2), nprimaries);
for ii = 1 : nprimaries
    % define function handle
    hf = @(x) x.data(controlMap(:, :, ii));

    % process
    outImg(:, :, ii) = blockproc(inImg(:, :, ii), pixels_per_dixel, hf);
end

end