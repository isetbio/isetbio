function y = sumOfTwoDblExponentials(coef, x)
% Compute sum of two double exponential functions, with no mean offset.
%
% Syntax:
%   y = sumOfTwoDblExponentials(coef, x)
%
% Description:
%    Compute sum of two double exponential functions, with no mean offset.
%       y = A1 * exp(-C1 * abs(x)) + A2 * exp(-C2 * abs(x));
%
% Inputs:
%    coef - A vector containing the required four coefficients in the order
%           described below:
%               1 - A1 - Variable to multiply the first exponent by
%               2 - C1 - The negative is multipled by the absolute value of
%                        x and then the exponent thereof is calculated
%               3 - A2 - Variable to multiply the second exponent by
%               4 - C2 - The negative is multiplied by the absolute value
%                        of x and then the exponent thereof is calculated
%    x    - The values on which to compute the double exponential
%
% Outputs:
%    y    - The calculated sum of the double exponential functions.
%
% Optional key/value pairs:
%    None.
%

% History:
%    01/08/16  dhb  Added commments
%              dhb  Changed name and what is computed to be more standard
%    12/04/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match Wiki.

% Do it
y = coef(1) .* exp(-coef(2) * abs(x)) + coef(3) .* exp(-coef(4) * abs(x));
