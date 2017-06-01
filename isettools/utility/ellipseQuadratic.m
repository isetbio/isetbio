function Q = ellipseQuadratic(ellipseParameters, rfRadius)
%ELLIPSEQUADRATIC
%
% Return the quadratic form for an ellipse based on the principal axes
% (y,x) and a rotation
%
%   bivariateGaussian =  exp(- (vec' * Q * vec)
% 
% It should be the case that when vec is a set of vectors that all satisfy
% 
%      1 = vec'*Q*vec
%
% they form the ellipsoid that is the 1 st dev ellipsoid for a bivariate
% Gaussian distribution.  You can find those by putting in all unit vectors
% and solving
%
%      v = unitVec'*Q*unitVec
%
% and then setting vec = unitVec/sqrt(v).  Then
%
%     vec' * Q * vec = unitVec'/sqrt(v) * Q unitVec/sqrt(v)
%                    = (1/v) * unitVec' * Q * unitVec
%                    = (1/v) * v = 1
%
% How to get unit vectors
%s
%  [~, pts] = ieShape('circle');
%
% 
% EJ lab code:
%   R = rmatrix2(angle / (2*pi) * 360);
%   L = params.sd_radius * [sd(1) 0; 0 sd(2)];
%   X = R * L * [x_circle; y_circle];
%         
% Example
%   eParams = [1 1 45];
%   Q = ellipseQuadratic(eParams);
%
% 
% JRG/BW 

D = diag(ellipseParameters(1:2));
R = eye(2);
% R = [cosd(ellipseParameters(3)) -sind(ellipseParameters(3));
%     sind(ellipseParameters(3))   cosd(ellipseParameters(3))];
Q = sqrt(D)*R;
% Q = R*(.5*rfRadius*D);        

end