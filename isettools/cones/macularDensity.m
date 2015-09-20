function density = macularDensity(eccDeg, varargin)
% Compute macular pigment optical density as a function of eccentricity
%    
%    density = macularDensity(eccentricity)
%
% Inputs:
%   eccDeg - eccentricity in degrees
%
% Outputs:
%   density - macular pigment optical density
%
% Notes:
%   1) Macular pigment density is roughly symmetric and thus we approximate
%      the 2D position by 1D eccentricity
%   2) The lorentzian function is fitted from data grabbed from figure 2(B)
%      in the reference paper. The data is stored in macularDensity.mat
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
if notDefined('eccDeg'), eccDeg = 0; end

% Compute density with the lorentz function
% Here, we force the model to be symmetric and have 0 density at infinite
% eccentricity
density = 1.261 ./ (eccDeg.^2 + 3.6028);

end