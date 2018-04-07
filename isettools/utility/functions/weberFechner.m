function y = weberFechner(coef, x)
% Calculate the Weber Fechner tvi function
%
% Syntax:
%   y = weberFechner(coef, x)
%
% Description:
%    Compute the Weber-Fechner sensitivity function used by Angueyra and
%    Rieke, 2013. to describe how rod flash response gain changes as a
%    function of background. This is their equation 1, with a log10
%    applied.
%       y = log10(1 ./ (1 + x / coef(1)));
%
%    They refer to this equation as a Weber-Fechner function, hence the
%    name of this function.
%
%    Examples are located within the code. To access the examples, type
%    'edit weberFechner.m' into the Command Window.
%    
% Inputs
%    coef - Coefficient of the equation, their I0.  This is the background
%           intensity that reduces the gain by a factor of 2 from its level
%           in the dark.
%    x    - Background intensities
%
%  Outputs:
%    y    - Calculated log10 of gain ratio gammaB/gammaD as a function of
%           background, where gammaB is the gain on a given background to
%           gain in darkness.  Gain is expressed as pA/R* where pA is the
%           photocurrent and R* is the number of isomerizations in the
%           flash.
% 
% Optional key/value pairs:
%    None.
%
% References:
%    Angueyra, Juan M., and Fred Rieke. "Origin and effect of
%    phototransduction noise in primate cone photoreceptors." Nature
%    neuroscience 16.11 (2013): 1692-1700.
%
% See Also:
%    v_osStepFlash
%

% History:
%    01/08/16  dhb  Added comment.
%    12/05/17  jnm  Formatting
%    01/06/18  dhb  Added comments that express our current lack of
%                   understanding.
%    01/26/18  jnm  Formatting update to match Wiki.

% Examples:
%{
   x = logspace(-1, 2, 50);
   y = weberFechner(1, x);
   vcNewGraphWin; loglog(x, y);
%}

% Do it.
y = log10(1 ./ (1 + x ./ coef(1)));

end