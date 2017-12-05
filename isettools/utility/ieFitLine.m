function [slope, offset] = ieFitLine(x, y, method)
% Fit a line to the data in x, y
%
% Syntax:
%   [slope, offset] = ieFitLine(x, y, [method])
%
% Description:
%    Solve least squares line y = a * x + b
%    Returns slope (a) and offset (b)
%
%    The data must be in the columns of x and y
%
% Inputs:
%    x      - Causation data
%    y      - Response data
%    method - (Optional) The method by which you fit the x and y data to a
%             line. Default is 'leastSquares'. The options are:
%               {'oneline', 'onelineleastsquares', 'leastsquares'}
%               {'multiplelines', 'multiplelinesleastsquares'}
%
% Outputs:
%    slope  - The slope of the line. 'a' in the formula above
%    offset - The vertical offset from the x axis. 'b' in the formula above
%
% Notes:
%    * This routine has not been used yet or debugged carefully.
%    * [Note: JNM - Aren't there existing built-in MATLAB functions that
%      accomplish this? Possibly as a part of statistics? - Optimization
%      toolbox contains multiple solutions to this.]
%    * [Note: XXX - TODO: There is a simpler formula. Must re-derive and
%      use it instead of this.]
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    11/30/17  jnm  Formatting, add example
%

% Examples:
%{
    % The results for the offset is a minutely non-zero value (-3.3307e-16)
    % as opposed to the zero slope for a direct line fit.
    [slo, offs] = ieFitLine(1:10,2:2:20)
%}

if notDefined('method'), method = 'leastSquares'; end
method = ieParamFormat(method);

nData = size(y, 2);
nSamples = size(y, 1);
   
% We can have one x-variable generating many y-variables, such as one
% exposure time list producing voltages at many pixels
if nData > 1 && size(x, 2) == 1
    % Could set up a flag and just store the x variable once, rather than
    % waste space like this
    x = repmat(x, 1, nData);
end

switch method
    case {'oneline', 'onelineleastsquares', 'leastsquares'}
        x = [x(:), ones(length(x), 1)];  y = y(:);
        val = pinv(x) * y;
        slope  = val(1);
        offset = val(2);
        % [Note: XXX - TODO: There is a simpler formula. Must re-derive and
        % use it instead of this.]
    case {'multiplelines', 'multiplelinesleastsquares'}

        % y = Ax, so we solve x = A\y;
        onesCol = ones(nSamples, 1);
        
        for ii=1:nData    
            thisX = [x(:, ii), onesCol]; 
            val = pinv(thisX) * thisY;
            slope(ii) = val(1);
            offset(ii) = val(2);
        end
        
    otherwise
        error('unknown method')
end

end