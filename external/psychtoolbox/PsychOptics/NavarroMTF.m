function mtf = NavarroMTF(s)
%NAVARROMTF  Compute the MTF measured by Navarror, Artal and Williams
%   mtf = NAVARRORMTF(s)
% 
%   Compute the MTF measured by Navarro et al.
%
%   Navarro, R., Artal, P., & Williams, D. R. (1993). 
%   Modulation transfer function of the human eye as a function of retinal 
%   eccentricity. Journal of the Optical Society of America A, 10, 201?212.
%
%   Spatial frequency passed in cycles/deg.
%
%   See also WILLIAMSMTF, OTFTOPSF, WILLIAMSRESTMTF, DIFFRACTIONMTF,
%   WILLIAMSTABULATEDPSF, PSYCHOPTICSTEST.

% 10/20/2918   npc		Wrote it.

mtf = 0.78*exp(-0.172*s) + 0.22*exp(-0.037*s);
