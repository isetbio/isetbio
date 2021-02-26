function density = eccDensity(obj, eccDeg, varargin)
% Compute macular pigment optical density as function of eccentricity
%
% Syntax:
%	density = eccDensity(eccDeg);
%
% Description:
%    Compute the macular pigment's peak optical density as a function of
%    its eccentricity.
%
%    This currently has foveal peak optical density hard coded at 0.35.
%
% Inputs:
%    eccDeg   - eccentricity in degrees.  May be a vector.
%
% Outputs:
%	 density  - macular pigment peak optical density, one entry for every
%	            entry of the input.
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

% History:
%   xx/xx/xx  npc  Wrote it.
%   12/15/20  dhb  Moved into its own file so we can get help on it.
%   12/20/20  npc  Made it an object method.

% Examples:
%{
    m = Macular();
    m.eccDensity(5);
%}

p = inputParser;
p.addParameter('eccDegs2', []);
p.parse(varargin{:});
eccDegs2 = p.Results.eccDegs2;

% Compute density with the lorentz function 
% 
% Here, we force the model to be symmetric and have 0 density
% at infinite eccentricity.
if (isempty(eccDegs2))
    density = obj.density * 3.6028 ./ (eccDeg .^ 2 + 3.6028);
else
    density = obj.density * 3.6028 ./ (eccDegs2 + 3.6028);
end



