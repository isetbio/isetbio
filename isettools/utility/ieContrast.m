function contrast = ieContrast(input,varargin)
% Convert a signal to a unit contrast
%
%    contrast = ieContrast(input,varargin)
%
% Input:
%   input - Usually a current signal input to a class
%
% Input parameters
%   maxC - max contrast, a number between 0 and 1
%
% Output:
%   contrast - The input converted to a contrast
%
% The signal mean is subtracted and the signal is scaled so that the entire
% range is 1 unit.  
%
% N.B. If the data are constant, then the return is all zero contrast.
%
% Example
%   ieContrast(input,'maxC',0.5)
%
% See also: computeSeparable.m
% 
% JRG/BW ISETBIO Team, 2016

%% Parse
p = inputParser;
p.addRequired('input',@isnumeric);
p.addParameter('maxC',1,@(x)((x >= 0) && ( x <= 1)));

p.parse(input,varargin{:});
input = p.Results.input;
maxC  = p.Results.maxC;

%% Find the range and mean, knock yourself out

range = max(input(:)) - min(input(:));
if range == 0
    warning('Constant data, hence zero contrast.');
    contrast = zeros(size(input));
else
    mn = mean(input(:));
    contrast = maxC*(input - mn)/range;
end

end
