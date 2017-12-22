function y = hill(coef, x)
% Compute hill function
%
% Syntax:
%   y = hill(coef, x)
%
% Description:
%    Compute the hill function
%        y = 1.0 ./ (1 + (coef(1) ./ x) .^ coef(2));
%
% Inputs:
%    coef - The coefficients outside of the ligand concentration required
%           to calculate the hill equation. They are:
%               1 - The ligand concentration producing half occupation
%               2 - The Hill coefficient, describing cooperativity (or
%                   possibly other biochemical properties, depending on the
%                   context in which the Hill equation is being used). The
%                   cooperativity is described below in the notes section.
%    x    - Free (unbound) ligand concentration
%
% Outputs:
%    y    - The fraction of the receptor protein concentration that is
%           bound to ligand.
%
% Notes:
%    * The hill coefficient's value describes the cooperativity of ligand
%      binding in the following way:
%         n > 1 - Positively cooperative binding: Once one ligand molecule
%                 is bound to the enzyme, its affinity for other ligand
%                 molecules increases. For example, the Hill coefficient of
%                 oxygen binding to haemoglobin (an example of positive
%                 cooperativity) falls within the range of 1.7-3.2.[2]
%         n < 1 - Negatively cooperative binding: Once one ligand molecule
%                 is bound to the enzyme, its affinity for other ligand
%                 molecules decreases.
%         n = 1 - Noncooperative (completely independent) binding: The
%                 affinity of the enzyme for a ligand molecule is not
%                 dependent on whether or not other ligand molecules are
%                 already bound. When n=1, we obtain a model that can be
%                 modeled by Michaelis–Menten kinetics,[7] in which
%                 K_D = K_A = K_M the Michaelis-Menten constant.

% History:
%    01/08/16  dhb  Added comment.
%    12/04/17  jnm  Formatting

% Do it
y = 1.0 ./ (1 + (coef(1) ./ x) .^ coef(2));
