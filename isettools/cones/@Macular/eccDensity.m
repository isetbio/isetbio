function density = eccDensity(eccDeg)
% Compute macular pigment optical density as function of eccentricity
%
% Syntax:
%	density = macularDensity(eccDeg);
%
% Description:
%    Compute the macular pigment's peak optical density as a function of
%    its eccentricity.
%
% Inputs:
%    eccDeg   - eccentricity in degrees.  May be a vector.
%
% Outputs:
%	 density  - macular pigment optical density, one entry for every entry
%	            of the input.
%
% Notes:
%   1) Macular pigment density is roughly symmetric and thus we
%      approximate the 2D position by 1D eccentricity
%   2) The lorentzian function is fitted from data grabbed from
%      figure 2(B) in the reference paper. The data is stored
%      in macularDensity.mat
%
% References:
%   Putnam, C. M., & Bland, P. J. (2014). Macular pigment
%   optical density spatial distribution measured in a subject
%   with oculocutaneous albinism. Journal of Optometry, 7(4),
%   241-245.
%

% Compute density with the lorentz function 
% 
% Here, we force the model to be symmetric and have 0 density
% at infinite eccentricity
density = 0.35 * 3.6028 ./ (eccDeg .^ 2 + 3.6028);

end
