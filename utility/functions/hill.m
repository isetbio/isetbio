function y = hill(coef, x)
% Compute hill function
%
% Syntax:
%     y = hill(coef, x)
%
% Description:
%    Compute the hill function
%        y = 1.0 ./ (1 + (coef(1) ./ x) .^ coef(2));
%
% Inputs:
%    coef - The coefficients outside of the ligand concentration required
%           to calculate the hill equation. They are:
%               1 - The ligand concentration K_A producing half occupation
%               2 - The Hill coefficient n, describing cooperativity (or
%                   possibly other biochemical properties, depending on the
%                   context in which the Hill equation is being used). The
%                   cooperativity is described below in the notes section.
%    x    - Free (unbound) ligand concentration
%
%    The Hill coefficient's value describes the cooperativity of ligand
%    binding in the following way:
%         n > 1 - Positively cooperative binding: Once one ligand molecule
%                 is bound to the enzyme, its affinity for other ligand
%                 molecules increases. For example, the Hill coefficient of
%                 oxygen binding to haemoglobin (an example of positive
%                 cooperativity) falls within the range of 1.7-3.2.
%                 [Weiss, 1997]
%         n < 1 - Negatively cooperative binding: Once one ligand molecule
%                 is bound to the enzyme, its affinity for other ligand
%                 molecules decreases.
%         n = 1 - Noncooperative (completely independent) binding: The
%                 affinity of the enzyme for a ligand molecule is not
%                 dependent on whether or not other ligand molecules are
%                 already bound. When n = 1, we obtain a model that can be
%                 modeled by Michaelis–Menten kinetics,[Alon, 2007] in
%                 which K_D = K_A = K_M the Michaelis-Menten constant.
% Outputs:
%    y    - The fraction of the receptor protein concentration that is
%           bound to ligand.
%
% Optional key/value pairs:
%    None.
%
% References:
%    https://en.wikipedia.org/wiki/Hill_equation_(biochemistry).  At least
%    some of the description above is taken verbatim from this site.
%
%    Weiss, J. N. (1 September 1997). "The Hill equation revisited: uses
%    and misuses". The FASEB Journal. 11 (11): 835?841.
%
%    Alon, Uri (2007). An Introduction to Systems Biology: Design
%    Principles of Biological Circuits ([Nachdr.] ed.). Boca Raton, FL:
%    Chapman & Hall.

% History:
%    01/08/16  dhb  Added comment.
%    12/04/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match Wiki.

% Do it
y = 1.0 ./ (1 + (coef(1) ./ x) .^ coef(2));
