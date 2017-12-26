function y = weberFechner(coef, x)
% Calculate the Weber Fechner tvi function
%
% Syntax:
%   y = weberFechner(coef, x)
%
% Description:
%    Compute the Weber-Fechner Threshold-versus-Intensity Step function
%       y = log10(1 ./ (1 + x / coef(1)));
%
% Inputs
%    coef - Coefficient
%    x    - Step intensities (scalar or vector)
%
%  Outputs:
%    y    - Calculated luminance
%
% Notes:
%    * [Note: BW - I do not understand this function.  I put in an
%    example and made it run.  What is this output?] 
%

% Example
%{
   x = logspace(-1,2,50);
   y = weberFechner(1,x);
   vcNewGraphWin; loglog(x,y);
%}
% History:
%    01/08/16  dhb  Added comment.
%    12/05/17  jnm  Formatting

% Do it.
y = log10(1 ./ (1 + x ./ coef(1)));

end