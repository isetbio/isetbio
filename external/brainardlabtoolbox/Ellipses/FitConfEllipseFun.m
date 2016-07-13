function [f]= FitConfEllipseFun(k,x,u,Q,percent)
% [f] = FitConfEllipseFun(k,x,u,Q,percent)
%
% Penalty function for finding confidence ellipse.
%
% 1/8/04    dhb     Wrote it.

usub = x-u(:,ones(1,size(x,2)));
values = EvaluateQuadratic(u,Q,x);
targetN = round(percent*size(x,2)/100);
actualN = length(find(values < k));
f = abs(targetN-actualN);