function contrast = ieContrast(input, varargin)
% Convert a signal to a contrast represenation
%
% Syntax:
%   contrast = ieContrast(input)
%
% Description:
%    Convert a signal to a a contrast representation. The signal mean is
%    subtracted and the signal is scaled so that the entire range is 1
%    unit.
%
%    N.B. If the data are constant, then the return is all zero contrast.
%
% Inputs:
%    input    - Usually a current signal input to a class
%
% Output:
%    contrast - The input converted to a contrast
%
% Optional key/value pairs:
%    'maxC'   - scalar between 0 and 1 (default 1).  Maximum
%               contrast of the output.  That is, contrast as defined above
%               is then scaled by this number.
%
% See Also:
%    computeSeparable.m
% 

% History:
%    xx/xx/16  JRG/BW  ISETBIO Team, 2016
%    11/30/17  jnm  Formatting

% Examples:
%{
    input = 1:50;
    ieContrast(input, 'maxC', 0.5)
%}

%% Parse
p = inputParser;
p.addRequired('input', @isnumeric);
p.addParameter('maxC', 1, @(x)((x >= 0) && ( x <= 1)));

p.parse(input, varargin{:});
input = p.Results.input;
maxC  = p.Results.maxC;

%% Find the range and mean, knock yourself out
range = max(input(:)) - min(input(:));
if range == 0
    warning('Constant data, hence zero contrast.');
    contrast = zeros(size(input));
else
    mn = mean(input(:));
    contrast = maxC * (input - mn) / range;
end

end
