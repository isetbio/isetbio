function [u,Q,k,encloseError] = FitConfEllipse(x,percent,u,K)
% [u,Q,k] = FitConfEllipse(x,percent,[u],[K])
%
% Find the confidence ellipse for the passed data.
% Data should be in columns of passed matrix x, where
% each column is an observation
%
% 1/8/04    dhb     Wrote it.
% 3/22/04   dhb     Allow optional passing of u, K.

% Compute mean vector
if (nargin < 3 | isempty(u))
    u = mean(x,2);
end
if (nargin < 4 | isempty(K))
    K = cov(x');
end

Q = inv(K);
nTest = 300;
testKs = linspace(1,10,nTest);
for i = 1:nTest
    errors(i) = FitConfEllipseFun(testKs(i),x,u,Q,percent);
end
[encloseError,minIndex] = min(errors);
k = testKs(minIndex);




