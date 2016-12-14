function [A,Ainv,Q,ellParamsFit] = EllipsoidFit(x,ellParams0,offset)
% [A,Ainv,Q,ellParamsFit] = EllipsoidFit(x,[ellParams0],[offset])
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

% Offset?
if (nargin < 3 || isempty(offset))
    OFFSET = false;
else
    OFFSET = true;
end

% Get infor on ellipse location
maxX = max(x,[],2);
minX = min(x,[],2);
meanX = mean(x,2);

% Have a go at reasonable initial values
if (nargin < 2 || isempty(ellParams0))
    ellRanges = maxX-minX;
    ellParams0 = [ellRanges' 0 0 0]';
end
vlb = [1e-3 1e-3 1e-3 0 0 0]';
vub = [1e3 1e3 1e3 2*pi 2*pi 2*pi]';

% Add on offset parameters if we're searching on an offset.
% Bound reasonably based on data range, to prevent really nutty
% fits.
if (OFFSET)

    
    % Use these to initialize and set bounds
    ellParams0 = [ellParams0 ; [meanX(1) meanX(2) meanX(3)]' ];
    vlb = [vlb ; [minX(1) minX(2) minX(3)]' ];
    vub = [vub ; [maxX(1) maxX(2) maxX(3)]' ];
end

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

% Handle offset case
if (length(ellParams) == 9)
    x(1,:) = x(1,:) - ellParams(7);
    x(2,:) = x(2,:) - ellParams(8);
    x(3,:) = x(3,:) - ellParams(9);
    ellParams = ellParams(1:6);
end

[A,Ainv,Q] = EllipsoidMatricesGenerate(ellParams);
vectorLengths = diag(x'*Q*x);
f = sqrt(mean((vectorLengths-1).^2));
end
