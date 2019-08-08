function R = deg2rotm(deg)
% Takes in an angle (degrees) and returns a 2x2 rotation matrix.
%
% Syntax:
%   R = deg2rotm(deg)
%
% Description:
%   Takes in an angle (degrees) and returns a 2x2 rotation matrix.
%
% Inputs:
%  deg              - The desired angle of rotation (in degrees)
%
% Outputs:
%   R               - A 2x2 rotation matrix 
%
% Optional key/value pairs:
%   None.
%
% Examples are provided in the source code.
%
% See also: rotm2deg, eul2rotm, rotm2eul
%

% History
%  8/23/18  mab  Created.

% Examples:
%{
    R = deg2rotm(33);
 
    % Check that we get an orthonormal matrix
    Icheck = R*R';
    I = eye(2);
    if (max(abs(Icheck-I)) > 1e-10)
        error('Did not produce an orthonornal matrix');
    end
%}

% Brute force construction.
R = [cosd(deg),-1*sind(deg);sind(deg),cosd(deg)];

end