function res = rgcExtractStdFromCovariance(covMatrix)
% Convert the covariance representation to a std dev. and angle
%
%   res = rgcExtractStdFromCovariance(covMatrix)
%
% res is a structure with
%
%  std:  In microns of the principal directions of the covariance ellipse
%  directions:  The vectors of the principal directions
%  angle: The angle of the principal direction.
%
%
% 2010 Vistasoft at Stanford

if ~rgcIsCovarianceMatrix(covMatrix)
    error('Covariance matrix is not 2x2 symmetric positive definite')
end

% Get the principal directions and eigenvalues of these directions from the
% covariance definition
[v,d] = eig(covMatrix);

% Convert to direction
res.std = [sqrt(d(1,1)) sqrt(d(2,2))];
res.directions = v;
res.angle = acos(v(1,1));

end