function density = macularDensity(ecc, varargin)
% Compute macular pigment optical density as a function of eccentricity
%    
%    density = macularDensity(eccentricity)
%
% Inputs:
%   ecc - eccentricity in degrees
%
% Outputs:
%   density - macular pigment optical density
%
% Notes:
%   1) Macular pigment density is roughly symmetric and thus we approximate
%      the 2D position by 1D eccentricity
%   2) The lorentzian function is fitted from data grabbed from figure 2(B)
%      in the reference paper
%
% Reference:
%   Putnam, C. M., & Bland, P. J. (2014). Macular pigment optical density
%   spatial distribution measured in a subject with oculocutaneous
%   albinism. Journal of Optometry, 7(4), 241?245.
%
% See also:
%   macularCreate, macularSet, macularGet
%
% HJ, ISETBIO TEAM, 2015

% Check inputs
if notDefined('ecc'), error('eccentricity required'); end

% Compute density with the lorentz function
density = 1.2594 ./ ((ecc - 0.0338).^2 + 3.5972);

end