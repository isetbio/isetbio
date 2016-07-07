function [A,Ainv,Q,ellParamsFit] = EllipsoidFit(x,ellParams0)
% [A,Ainv,Q,ellParamsFit] = EllipsoidFit(x,[ellParams0])
%
% Find the ellipsoid that goes through a set of passed points in the
% columns of matrix x.  The method follows that described in:
%   Poirson AB, Wandell BA, Varner DC, Brainard DH. 1990. Surface
%   characterizations of color thresholds. J. Opt. Soc. Am. A 7: 783-89
% 
% You can pass an initial set of parameters.  If this argument is missing
% or empty, the routine tries to do something sensible.
% 
% See EllipsoidMatricesGenerate for the parameter convention.
% 
% 7/4/16  dhb  Wrote it.



% Set reasonable bounds on parameters
% Have a go at reasonable initial values
if (nargin < 2 || isempty(ellParams0))
    ellRanges = max(x,[],2)-min(x,[],2);
    ellParams0 = [ellRanges' 0 0 0]';
end
vlb = [1e-3 1e-3 1e-3 0 0 0]';
vub = [1e3 1e3 1e3 2*pi 2*pi 2*pi]';

%% Fit that sucker
%
% I coded up the global search method, but it is very slow compared with
% fmincon alone, and fmincon seems to be fine.
USEGLOBAL = false;
if (~USEGLOBAL)
    options = optimset('fmincon');
    options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off','Algorithm','active-set');
    ellParamsFit = fmincon(@(ellParams)FitEllipseFunction(ellParams,x),ellParams0,[],[],[],[],vlb,vub,[],options);
else
    opts = optimoptions(@fmincon,'Algorithm','interior-point');
    problem = createOptimProblem('fmincon','objective', ...
        @(ellParams)FitEllipseFunction(ellParams,x),...
        'x0',ellParams0,'lb',vlb,'ub',vub,'options',opts);
    gs = GlobalSearch;
    ellParamsFit = run(gs,problem);
end

[A,Ainv,Q] = EllipsoidMatricesGenerate(ellParamsFit);


end

%% Error function for the fit
function f = FitEllipseFunction(ellParams,x)

[A,Ainv,Q] = EllipsoidMatricesGenerate(ellParams);
vectorLengths = diag(x'*Q*x);
f = sqrt(mean((vectorLengths-1).^2));
end
