function [im, mn, mx] = ieScale(im, b1, b2)
% Scale the value in im into a specified range
%
% Syntax:
%   [im, mn, mx] = ieScale(im, [b1], [b2])
%
% Description:
%    Scale the values in im into the specified range. There can be one or
%    two bounds provided. Changes in the provided syntax will produce
%    different scaling operations
%
%   You will scale the values in im into the specified range. There can be
%   one or two bounds. The calculation is performed by scaling the data
%   between 0 and 1 and then applying b1 and b2 if they are provided.
%   	(im - min) / (max - min), and then (b2 - b1) * im + b1;
%
%    Examples contained in the code.
%
% Inputs:
%    im - Input data you wish to scale. If no bounds are provided, data
%         will scale between 0 and 1
%    b1 - (Optional) If only argument, the upper bound, else the lower
%    b2 - (Optional) If provided, this is the upper bound for scaling
%
% Outputs:
%    im - The scaled data
%    mn - The original lower bound
%    mx - The original upper bound
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/22/17  jnm  Formatting
%    01/16/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    im = -10:50
    im = ieScale(-10:50)
%}
%{
    im = -10:50
    im = ieScale(-10:50, 25)
%}
%{
  im = -10:50;
  [min(im(:)), max(im(:))]
  im = ieScale(im, 20, 90);
  [min(im(:)), max(im(:))]
%}

if notDefined('im'), error('Input data must be defined.'); end

% Find data range
mx = max(im(:));
mn = min(im(:));

% If only one bounds argument, just set peak value
if nargin == 2
    im = im * (b1 / mx);
    return;
end

% If 0 or 2 bounds arguments, we first scale data to 0, 1
im = (im - mn) / (mx - mn);

if nargin == 1
    % No bounds arguments, assume 0, 1
    b1 = 0;
    b2 = 1;
elseif nargin == 3
    if b1 >= b2
        error('ieScale: bad bounds values.');
    end
end

% Put the (0, 1) data into the range
range = b2 - b1;
im = range * im + b1;

end