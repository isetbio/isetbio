function r = GetMahalanobisDistance(x,u,K)
% r = GetMahalanobisDistance(x,u,K)
% 
% Compute the MahalanobisDistance for the passed
% columns of x.  Formula from Duda and Hart, p. 24.
%
% 6/30/05   dhb     Wrote it.

n = size(x,2);
r = zeros(1,n);
Q = inv(K);
for i = 1:n
    r(i) = sqrt((x(i)-u)'*Q*(x(i)-u));
end
