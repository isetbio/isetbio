function cMatrx = rgcCovarianceMatrix(sxx, syy, theta)
% Builds a covariance matrix with sxx and syy as specified, then rotates it
% with an angle theta (in radians)
%
% Typically we are setting up values in cone samples, not microns.
%
% Use sdInMicrons/pixelGet to convert from sd in microns to sd in cone
% counts.
%
% Theta is presumably in radians.
%
% Sd = rgcCovarianceMatrix(sxx,syy,theta);
%
% Ex:
%  Sd = rgcCovarianceMatrix(1, 10, pi/8);
%  support = [25,25]; spacing = [1,1];
%  g = getGaussian(support, spacing, Sd);imagesc(g);

if notDefined('sxx'), error('Std Dev in x needed');  end
if notDefined('syy'), error('Std Dev in y needed');  end
if notDefined('theta'), error('Rotation angle needed');  end

% Std Dev matrix in (x,y)
Sd = [sxx 0; 0 syy];

% Rotation matrix
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Build covariance matrix using rotation and standard deviations
cMatrx = (Sd*R)'*(Sd*R);

end