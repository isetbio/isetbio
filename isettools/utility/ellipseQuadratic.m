function Q = ellipseQuadratic(ellipseParameters, rfRadius)
% Deprecated.  This seems wrong to BW in many ways.  And not used elsewhere
% So, removing.
%
% OLD: Quadratic form of an ellipse based on principal axes (y, x) & a rotation
%
% Syntax:
%   Q = ellipseQuadratic(ellipseParameters, rfRadius)
%
% Description:
%    Return the quadratic form, Q, for an ellipse based on the principal axes
%    (y,x) and a rotation
%
%           bivariateGaussian =  exp(-(vec' * Q * vec))
%
%    It should be the case that when vec is a set of vectors that all
%    satisfy the formula
%
%           1 = vec' * Q * vec
%
%    They form the ellipsoid that is the 1 st dev ellipsoid for a bivariate
%    Gaussian distribution.  You can find those by putting in all unit
%    vectors and solving
%
%           v = unitVec' * Q * unitVec
%
%    and then setting vec = unitVec/sqrt(v). Then
%
%           vec' * Q * vec = unitVec' / sqrt(v) * Q unitVec / sqrt(v)
%                          = (1 / v) * unitVec' * Q * unitVec
%                          = (1 / v) * v = 1
%
%    How to get unit vectors
%
%           [~, pts] = ieShape('circle');
%
% Inputs:
%    ellipseParameters - The ellipse parameters
%    rfRadius          - (Optional) Radius
%
% Outputs:
%    Q                 - The quadratic form of the ellipse
%
% Notes:
%    * [Note: JNM - rfRadius is never used. Why is it still included?]
%    * [Note: JNM - This file itself is the only reference to this function
%      that I can find. Is there a particular reason for keeping it?]
%

% History:
%    xx/xx/xx  JRG/BW  Created
%    12/01/17  JNM     Formatted

% Examples:
%{
    eParams = [1 1 45];
	Q = ellipseQuadratic(eParams);
%}

% EJ lab code:
%    R = rmatrix2(angle / (2 * pi) * 360);
%    L = params.sd_radius * [sd(1) 0; 0 sd(2)];
%    X = R * L * [x_circle; y_circle];
error('deprecated');
end
%{
% Original code
D = diag(ellipseParameters(1:2));
R = eye(2);
% R = [cosd(ellipseParameters(3)) -sind(ellipseParameters(3));
%      sind(ellipseParameters(3))  cosd(ellipseParameters(3))];
Q = sqrt(D) * R;
% Q = R * (.5 * rfRadius * D);        
%}
end