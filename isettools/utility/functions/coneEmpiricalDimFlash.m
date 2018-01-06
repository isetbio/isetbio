function response = coneEmpiricalDimFlash(coef, t)
% Return estimate of the Schnapf et al. dim flash cone photocurrent response.
%
% Syntax:
%   response = coneEmpiricalDimFlash(coeffs, t)
%
% Description:
%    Using the coefficients and time provided by the user, calculate an
%    emprically-derived description of the dim flash cone photocurrent
%    response, following the formula in the Schnapf et al., 1990.
%
%    The fit uses a damped oscillator with S-shaped onset.
%
%    The calculation is as follows:
%       fit = scFact * (((t / TR) ^ 3) / (1 + ((t / TR) ^ 3))) ...
%             * exp(-((t / TD))) ...
%             * cos(((2 * pi * t) / TP) + (2 * pi * Phi / 360));
%
%    See the papter, Table 2 (p. 693) for some coefficient values.
%
% Inputs:
%    coef     - The coefficients for the calculations. They will follow the
%               following order in the provided vector:
%                   1 - scFact - Scaling Factor
%                   2 - TauR   - Rising Phase Time Constant
%                   3 - TauD   - Damping Time Constant
%                   4 - TauP   - Oscillator Period
%                   5 - Phi    - Oscillator Phase
%    t        - Time (in milliseconds)
%
% Outputs:
%    response - The calculated flash response
%
% References:
%    * Schnapf et al., 1990, Visual transduction in cones of the monkey
%      Macaca Fascicularis, p. 693.  Paper available at
%      https://goo.gl/EfqVwK

% History:
%    04/xx/11  Angueyra  Created 
%    04/xx/11  Rieke     Replaced gaussian by exponential decay
%    01/08/16  dhb       Rename, clean
%    12/04/17  jnm       Formatting

%% Give names to the coefficient vector entries
ScFact = coef(1);  % Scaling Factor
TauR = coef(2);    % Rising Phase Time Constant
TauD = coef(3);    % Dampening Time Constant
TauP = coef(4);    % Oscillator Period
Phi = coef(5);     % Oscillator Phase

%% Compute the response
response = ScFact .* (((t ./ TauR) .^ 3) ./ (1 + ((t ./ TauR) .^ 3))) ...
    .* exp(-((t ./ TauD))) ...
    .* cos(((2 .* pi .* t) ./ TauP) + (2 * pi * Phi / 360));

%% Some earlier versions now not used
% response = ScFact .* (((t ./ TauR) .^ 3) ./ (1 + ((t ./ TauR) .^ 3))) ...
%       .* exp(-((t ./ TauD) .^ 2)) ...
%       .* cos(((2 .* pi .* t) ./ TauP) + (2 * pi * Phi / 360));
% response = ScFact .* (((t ./ TauR) .^ 4) ./ (1 + ((t ./ TauR) .^ 4))) ...
%       .* exp(-((t ./ TauD))) ...
%       .* cos(((2 .* pi .* t) ./ TauP) + (2 * pi * Phi / 360));
