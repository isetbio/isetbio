function gaussianDistribution =  getGaussian(rfSupport, rfCov)
% Creates a 2D gaussian distribution
%
% Syntax:
%	gaussianDistribution =  getGaussian(rfSupport, rfCov)
%
% Description:
%    This routine generalizes the fspecial call because the bivariate
%    Gaussian can have a general layerariance matrix.
%
% Inputs:
%    rfSupport            - The support structure, containing x and y.
%    rfCov                - A 2 x 2 layer covariance matrix
%
% Outputs:
%    gaussianDistribution - The 2D Gaussian distribution
%
% cov is a 2x2 layer covariance matrix.
%
% Notes:
%    * [Note: the meshgrid in the example below does not require both
%      -10:10 inputs to create the desired output.]

% History:
%    xx/xx/09       (c) Stanford Synapse Team 2009
%    12/15/17  jnm  Formatting

% Examples:
%{
    [rfS.X rfS.Y] = meshgrid(-10:10, -10:10);
    rfCcov = [5 0; 0 9]; 
    g = getGaussian(rfS, rfCcov);
    vcNewGraphWin;
    mesh(rfS.X, rfS.Y, g);
%}

%% N-dimensional Gaussian
% We use this instead of fspecial because we want to allow asymmetric
% covariance matrices (non-rotational)
% M is crazy, too many coordinates
% S is the inverse of the covariance matrix.
gaussFunc = @(M, S)(1 / (2 * pi * sqrt(det(rfCov)))) ...
    * exp(-(1 / 2) * (M') * (rfCov \ M));

%% Check input
if notDefined('rfSupport'), error('rfSupport.X / .Y required'); end
if notDefined('rfCov'), error('rfCov required'); end

% Make sure the covariance matrix has sensible properties
if any(size(rfCov) ~= [2 2]), error('rfCov has to be a 2x2 matrix'); end
if (det(rfCov) == 0), error('rfCov is not invertible'); end
if (trace(rfCov) <= 0 || det(rfCov) < 0)
    error('rfCov is not positive');
end
if (abs(rfCov(1, 2) - rfCov(2, 1)) > 1e-6)
    error('rfCov is not symmetric');
end

%% Perform calculation

% Apply the gaussian functions
M = [rfSupport.X(:) rfSupport.Y(:)]';
gaussianDistribution = gaussFunc(M, rfCov); % See inline function above

% We wasted a lot of time computing the off-diagonal entries of(M inv(S) M)
gaussianDistribution = diag(gaussianDistribution);

% reshape and making sure the volume is 1
gaussianDistribution = reshape(gaussianDistribution, ...
    size(rfSupport.X, 2), size(rfSupport.Y, 1));
gaussianDistribution = gaussianDistribution / sum(gaussianDistribution(:));
% vcNewGraphWin;
% mesh(rfSupport.X, rfSupport.Y, gaussianDistribution)

end