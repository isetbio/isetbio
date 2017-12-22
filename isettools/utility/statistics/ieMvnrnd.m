function s = ieMvnrnd(mu, Sigma, K)
% Multivariate normal random number generator
%
% Syntax:
%   s = ieMvnrnd(mu, Sigma, K)
%
% Description:
%    This is a multivariate normal random number generator.
%
%    The data will be d dimensional.  n is the number of samples, although
%    this can be specified two ways (see below).
%
% Inputs:
%    mu    - (Optional) A matrix of n x d, that is, n samples, by d
%            dimensions. Default is the scalar value 0.
%    Sigma - (Optional) A covariance matrix of d x d. Default is the scalar
%            value 1.
%    K     - The total number of samples. Typically used when mu contains
%            one row, at which point there is a repmat which rewrites mu to
%            K rows for K samples.  Matlab calls K the number of 'cases'.
%
% Outputs:
%    s     - The requested random numbers
%
% Notes:
%    * [Note: XXX - Draw n random d-dimensional vectors from a multivariate
%      Gaussian distribution with mean mu (n x d) and covariance matrix
%      Sigma (d x d).  K specifies how many samples for each condition.]
%    * [Note: JNM - Do we want to include a check for if K is defined?]
%

% History:
%    05/xx/04  XXX?   Add Third arg.
%    xx/xx/14  HJ/BW  ISETBIO TEAM, 2014
%    12/13/17  JNM    Formatting
% Examples:
%{
    % Univariate
    mu = zeros(1000, 1);
    Sigma = 1;
    s = ieMvnrnd(mu, Sigma);
    vcNewGraphWin;
    hist(s, 50);
%}
%{
    % Univariate
    mu = 0;
    Sigma = 1;
    K = 1000;
    s = ieMvnrnd(mu, Sigma, K); 
    vcNewGraphWin;
    hist(s, 50);
%}
%{
    % Bivariate Covariance
    mu = zeros(500, 2);
    Sigma = [.5 -.3;-.3 1];
    s = ieMvnrnd(mu, Sigma); 
    vcNewGraphWin;
    plot(s(:, 1), s(:, 2), '.'); 
    axis equal;
    set(gca, 'xlim', [-5 5], 'ylim', [-5 5])
%}

%% Argument
if notDefined('mu'), mu = 0; end
if notDefined('Sigma'), Sigma = 1; end

% If mu is column vector and Sigma not a scalar then assume user didn't
% read help but let them off and flip mu.
if size(mu, 2) == 1 && ~isscalar(Sigma), mu = mu'; end

% May 2004 take a third arg, cases. Makes it more compatible with Matlab's.
if nargin == 3, mu = repmat(mu, K, 1); end
[n, d] = size(mu);

if any(size(Sigma) ~= [d, d])
    error('Sigma must have dimensions d x d where mu is n x d.');
end

% Check for stats toolbox
if checkToolbox('Statistics Toolbox')
    % Matlab toolbox version is present. Use it.
    s = mvnrnd(mu, Sigma);
    return
end

try
    U = chol(Sigma);
catch
    [E, Lambda] = eig(Sigma);
    if (min(diag(Lambda)) < 0)
        error('Sigma must be positive semi-definite.');
    end
    U = sqrt(Lambda) * E';
end

s = randn(n, d) * U + mu;

end

%%
% {{{ END OF CODE --- Guess I should provide an explanation:
% 
% We can draw from axis aligned unit Gaussians with randn(d)
%     x ~ A * exp(-0.5 * x' * x)
% We can then rotate this distribution using
%     y = U' * x
% Note that
%     x = inv(U') * y
% Our new variable y is distributed according to:
%     y ~ B * exp(-0.5 * y' * inv(U' * U) * y)
% or
%     y ~ N(0, Sigma)
% where
%     Sigma = U' * U
% For a given Sigma we can use the chol function to find the corresponding
% U, draw x and find y. We can then adjust for a non-zero mean by just
% adding it on.
% 
% But the Cholsky decomposition function doesn't always work...
% Consider Sigma = [1 1;1 1]. Now inv(Sigma) doesn't actually exist, but
% Matlab's mvnrnd provides samples with this covariance st x(1)~N(0, 1)
% x(2) = x(1). The fast way to deal with this would do something similar to
% chol but be clever when the rows aren't linearly independent. However, I
% can't be bothered, so another way of doing the decomposition is by
% diagonalising Sigma (which is slower but works).
% if
%     [E, Lambda] = eig(Sigma)
% then
%     Sigma = E * Lambda * E'
% so
%     U = sqrt(Lambda) * E'
% If any Lambdas are negative then Sigma just isn't even positive
% semi-definite so we can give up.
% }}}