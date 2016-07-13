function [vals] = EvaluateQuadratic(u,Q,x)
% [vals] = EvaluateQuadratic(u,Q,x);
%
% Evaluate the quadratic form (x-u)'*Q*(x-u).
%
% 1/8/04    dhb     Wrote it.

usub = x-u(:,ones(1,size(x,2)));
vals = diag( usub'*Q*usub )';
