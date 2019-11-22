function [img, params] = image1d(pattern, varargin)
% Make a 2D image from a 1D pattern.
%
% Syntax:
%   img = image1d(pattern, [varargin])
%
% Description:
%    Use a 1D pattern to create a 2D image.
%
%    This function contains examples of usage. To access, simply enter
%    'edit image1d.m' into the Command Window.
%
% Inputs:
%    pattern - Matrix. A matrix representing the 1D pattern.
%
% Outputs:
%    img     - Matrix. A matrix representation of the image.
%    params  - Struct. A structure containing the image paramters.
%
% Optional key/value pairs:
%    rows    - Numeric. The row length.
%    rgb     - Matrix. A 1x3 matrix representing the RGB colors. This is
%              designed to leave the mean unchanged.
%    mn      - Numeric. Treat this value as the image mean.
%    pad     - Pad with mean, or maybe it will be edge?
%              (This is not yet implemented)
%

% Examples:
%{
    pattern = ones(1, 129) * 0.5;
    pattern(65) = 1;
    img = image1d(pattern); 
    vcNewGraphWin;
    imshow(img)

    x = (-63:64) / 128;
    f = 2;
    pattern = 0.5 * cos(2 * pi * f * x) + 0.5;
    [img, p] = image1d(pattern, 'rows', 256);
    imshow(img);

    [img, p] = image1d(pattern, 'rgb', [.2 .2 1], 'rows', 64);
    imshow(img);
%}

% Is there a matlab parameter parser we should be using?

% Init parameters
params.rgb = [1 1 1];
params.rows = length(pattern);
params.pad = [0 0];   % Extend the pattern

for ii = 1:2:length(varargin)
    thisP = ieParamFormat(varargin{ii});
    switch thisP
        case 'rgb'
            params.rgb = varargin{ii + 1};
        case 'rows'
            params.rows = varargin{ii + 1};
        case 'mean'
            params.mean = varargin{ii + 1};
        otherwise
            error('Unknown parameter %s\n', thisP);
    end
end

%% Make the image
% 1D image
img = repmat(pattern(:)', [params.rows, 1]);

% Extend the pattern here?
% Look at params.pad

% Color the contrast in RGB space, leaving the mean unchanged.

% Perhaps we should color the ccontrast, and not the mean?
% We could pull out the mean, multiply the contrast, and then add them
% together.
[img, r, c] = RGB2XWFormat(img);
if ~isfield(params, 'mean'), mn = mean(img); else, mn = params.mean; end
contrast = img - mn;
contrast = contrast * params.rgb(:)';
img = mn + contrast;
img = XW2RGBFormat(img, r, c);

end
