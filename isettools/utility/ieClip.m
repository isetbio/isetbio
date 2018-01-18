function im = ieClip(im, lowerBound, upperBound)
% Clip data to range specified arguments.
%
% Syntax:
%   im = ieClip(im, [lowerBound], [upperBound])
%
% Description:
%    Various types of clipping of data are supported. These include the
%    following formats.
%       * Providing a single variable - the data - data returned between
%         the range of 0 to 1
%       * Providing two variables - the data and a bound - sets the
%         boundary from the negative version of the bound to the positive
%         version of the bound
%       * Providing three variables - the data, an upper, and a lower
%         bound - The data will be bounded at the bottom by the lower bound
%         (2nd variable), and at the top by the upper bound (3rd variable),
%         there are some special situations:
%           - If the second variable is [] - no lower bound ('blank' var 2)
%           - If the third variable is [] - no upper bound ('blank' var 3)
%           - If both are [] (blank) - no bounding will occur
%       * Providing the data and a lower and an upper bound
%
%    Examples are included within the code.
%
% Inputs:
%    im         - The data you wish to clip
%    lowerBound - (Optional) If upperBound exists, the lowerBound, else the
%                 +/- boundary value for the data. Default 0.
%    upperBound - (Optional) The upper bound for the data. Default 1.
%
% Outputs:
%    im         - The clipped data
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    11/30/17  jnm  Formatting
%    12/21/17  BW   Fixed note by using fprintf and fixed exist()
%    01/16/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    im = 2 * randn([5, 5])
    ieClip(im, [], 1)
    ieClip(im, 0, [])
    ieClip(im, 1.35789)
%}
%{
    im = 10 * randn([5, 5])
    bound = 7.5;
    ieClip(im, [], 5) % sets the upper bound to 5, no lower bound
    ieClip(im, 0, 1)    % sets the lower to 0 and upper to 1
    ieClip(im, 0, [])   % sets the lower to 0, no upper bound
    ieClip(im)          % defaults to 0 1 range
    ieClip(im, bound)   % sets bound to +/- bound
%}

% Set up variables
if nargin == 1
    % Only im sent in. Default is [0, 1]
    lowerBound = 0;
    upperBound = 1;
    fprintf('ieClip:  Setting range to 0 1');
elseif nargin == 2
    % Reads this as [-l, l]
    lowerBound = -abs(lowerBound);
    upperBound = abs(lowerBound);
    fprintf('ieClip:  Setting range to [%.3e, %.3e]', lowerBound, ...
        upperBound);
end

if ~(~exist('lowerBound', 'var') || isempty(lowerBound))
    im(im<lowerBound) = lowerBound;
end

if ~(~exist('upperBound','var') || isempty(upperBound))
    im(im > upperBound) = upperBound;
end

end
