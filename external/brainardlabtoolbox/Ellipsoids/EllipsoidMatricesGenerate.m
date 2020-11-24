function [A,Ainv,Q] = EllipsoidMatricesGenerate(ellParams,varargin)
% Generate matrix description of ellipsoid from parameter vector
% 
% Syntax:
%     A = EllipsoidMatricesGenerate(ellParams)
% 
% Description:
%     Generate the 3 by 3 matrix A from a 6-vector of parameters.  When applied to points on 
%     the ellipsoid it specifies, the matrix A maps these points to the unit sphere in 3D.
%     Thus the matrix inv(A) maps points on the unit sphere to points on the ellipsoid.
%
%     Can also generate 2 by 2 matrices from a 3-vector of parameters.
% 
%     The column vector ellParams are the the diagonal entries of a 3 by 3 matrix D
%     and the Euler rotation angles (in radians) for a 3D rotation matrix
%     R.  In the 2D case, there are just two diagnonal entries and one
%     rotataion angle.
% 
%     The parameter matrix can optionally by 9 entries long, in which case the
%     last three entries are the coordinates of the center of the ellipsoid.
% 
%     The matrix D stretches the x, y, and z axes in the coordinate system of the unit sphere,
%     producing the three princple axes of the ellipsoid aligned to x, y,
%     z.  The larger the entry of D, the smaller the axis.
% 
%     The rotation matrix R rotates these axes to their desired orientations.
% 
%     Q is given as A'*A, and has the property that for points x on the ellipsoid,
%     x'*Q*x = k.  The constant determines the scale of the
%     ellipsoid/ellipse.
% 
%     The Euler angles are passed to eul2rotm and interpretted in its default
%     'ZYX' order.  Thus the Euler angles are in radians.
% 
%     The parameterization of Q follows that in 
%       Poirson AB, Wandell BA, Varner DC, Brainard DH. 1990. Surface
%       characterizations of color thresholds. J. Opt. Soc. Am. A 7: 783-89.
%     See particularly pp. 784-785.
%
%     We had some trouble deciding the right convention for angles.  For
%     the two dimensional case, this is set up so that the major axis of
%     the ellipse is rotated counter-clockwise relative to the x-axis
%     according to the passed angle.  Need to check more carefully, and
%     particularly for the 3-D case.
%
% Optional key/value pairs:
%    'dimension'    - What dimension is the ellipsoid? Can be 2 or 3.
%                     Default is 3. Passing 2 means an ellipse.

% History:
%   06/27/16  dhb  Back to the future.  Wrote this.  It feels like 1988.
%   08/16/18  dhb  Change parameterization to match paper.
%   11/20/18  dhb  Remove transpose from scalar arg to deg2rotm.
%   11/24/18  dhb  Flip convention of V, V' to match intuition about
%                  angles.  Checked for two-dimensional case.
%   11/19/20  dhb  Remove transpose from definition of V.  We now think
%                  that should not have been there.

% Examples:
%{
  ellParams = [1 0.5 45];
  [A,Ainv,Q] = EllipsoidMatricesGenerate(ellParams,'dimension',2);
%}

% Parse input
p = inputParser;
p.addRequired('ellParams',@isnumeric);
p.addParameter('dimension',3,@(x) (isnumeric(x) & isscalar(x)));
p.parse(ellParams,varargin{:});

% Convert between parameter vector and the S and V matrices
%
% If an offset was passed, we ignore it because we are just
% fitting 
switch(p.Results.dimension)
    case 2
        S = diag(ellParams(1:2));
        V = deg2rotm(ellParams(3));  
        if (length(ellParams) == 5)
            error('Translation not yet implemented');
        end
    case 3
        S = diag(ellParams(1:3));
        V = eul2rotm(deg2rad(ellParams(4:6)'));  
        if (length(ellParams) == 9)
            error('Translation not yet implemented')
        end
    otherwise
        error('Can only deal with dimension set to 2 or 3');
end

% Get A, Ainv and Q
A = S*V';
Ainv = inv(A);
Q = A'*A;

