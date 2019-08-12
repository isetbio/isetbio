function deg = rotm2deg(R)
% Extract rotation in degrees from 2D rotation matrix
%
% Syntax:
%   deg = rotm2deg(R)
%
% Description:
%   Takes in a 2x2 rotation matrix and returns an angle in degrees.
%
% Inputs:
%  R               - A 2x2 rotation matrix
%
% Outputs:
%  deg             - The desired angle of rotation (in degrees) 
%                    Returned angles wrap to [0,360].
%
% Optional key/value pairs:
%   None.
%
% Examples are provided in the source code.
%
% See also: deg2rotm, eul2rotm, rotm2eul
%

% History
%  8/23/18  mab  Created.

% Examples:
%{
    degs = [0 22 90 94 180 187 270 275];
    for ii = 1:length(degs)
       R = deg2rotm(degs(ii));
       deg1 = rotm2deg(R);
        if (abs(degs(ii)-deg1) > 1e-6)
            error('Failure to self-invert');
        end
    end
%}

% Handke the four quadrants by brute force.
% Angles come back in range [0,360];
cosTheta = R(1,1);
sinTheta = R(2,1);

if (cosTheta >= 0 & sinTheta >= 0)
    deg = acosd(cosTheta);
elseif (cosTheta >= 0 & sinTheta < 0)
    deg = asind(sinTheta)+360;
elseif (cosTheta < 0 & sinTheta >= 0)
    deg = acosd(cosTheta);
elseif (cosTheta < 0 & sinTheta < 0)
    deg = acosd(cosTheta)+ 2*(180-acosd(cosTheta));
else
    error('My computer is acting strange');
end

% Check
cosThetaCheck = cosd(deg);
sinThetaCheck = sind(deg);
if (abs(cosThetaCheck-cosTheta) > 1e-6 | abs(sinThetaCheck-sinTheta) > 1e-6)
    error('Did not extract angle properly. Fix the code.');
end

% If we were passed a rotation matrix, 
Rcheck = deg2rotm(deg);
if (abs(Rcheck(:)-R(:)) > 1e-6)
    error('Passed matrix does not appear to be a rotation');
end

end