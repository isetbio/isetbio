function [A,Ainv,Q,ellParamsFit] = EllipsoidFit(x,ellParams0,fitCenterOffset,isXYEllipse)
% Fit an ellipsoid (or ellipse) to data
%
% Syntax
%     [A,Ainv,Q,ellParamsFit] = EllipsoidFit(x,[ellParams0],[offset])
%
% Description:
%     Find the ellipsoid that goes through a set of passed points in the
%     columns of matrix x.  The method follows that described in:
%        Poirson AB, Wandell BA, Varner DC, Brainard DH. 1990. Surface
%        characterizations of color thresholds. J. Opt. Soc. Am. A 7: 783-89
% 
%     You can pass an initial set of ellipse parameters.  If this argument is missing
%     or empty, the routine tries to do something sensible.
%
%     If isXYEllipse is true, the z dimension center position is locked to
%     zero.  This is useful when we are using the ellipsoid code to fit data
%     which is restricted to the xy plane.  This is a bit of hack, and may not
%     be completely robust.
%
% Inputs:
%
% Outputs:
%
% See EllipsoidMatricesGenerate for the parameter convention.
% 
% 07/04/16  dhb  Wrote it.
% 08/16/18  dhb  Update initial guess to match new ellipsoid parameterization
%                Also allow angles to run from -2*pi to 2*pi, so search can
%                wrap around if it needs to.
% 04/22/19  dhb  Respect passed fitCenterOffset - this was being forced to
%                true when passed.

% Offset?
if (nargin < 3 || isempty(fitCenterOffset))
    fitCenterOffset = false;
end

if (nargin < 4 || isempty(isXYEllipse))
    isXYEllipse = false;
end

% Get info on ellipse location
maxX = max(x,[],2);
minX = min(x,[],2);
meanX = mean([maxX minX],2);
meanX = mean(x,2);

% Have a go at reasonable initial values
if (nargin < 2 || isempty(ellParams0))
    ellRanges = maxX-minX;
    ellParams0 = [1./ellRanges.^0.5' 0 0 0]';
    ellParams0(isinf(ellParams0)) = 1;
end
vlb = [1e-3 1e-3 1e-3 -2*pi -2*pi -2*pi]';
vub = [1e3 1e3 1e3 2*pi 2*pi 2*pi]';

% Add on offset parameters if we're searching on an offset.
% Bound reasonably based on data range, to prevent really nutty
% fits.
if (fitCenterOffset)
    
    % Set initial values and bounds on center
    rectFactor = 2;
    minFactor = 1;
    maxFactor = rectFactor-1;
    meanFactor = 0.05;
    if (~isXYEllipse)
        ellParams0 = [ellParams0 ; [meanX(1) meanX(2) meanX(3)]' ]; 
        vlb = [vlb ; [minX(1)+minFactor*(maxX(1)-minX(1))/rectFactor minX(2)+minFactor*(maxX(2)-minX(2))/rectFactor minX(3)+minFactor*(maxX(3)-minX(3))/rectFactor]' ];
        vub = [vub ; [minX(1)+maxFactor*(maxX(1)-minX(1))/rectFactor minX(2)+maxFactor*(maxX(2)-minX(2))/rectFactor minX(3)+maxFactor*(maxX(3)-minX(3))/rectFactor]' ];
    else
        ellParams0 = [ellParams0 ; [meanX(1) meanX(2) 0]' ];
        %vlb = [vlb ; [minX(1)+minFactor*(maxX(1)-minX(1))/rectFactor minX(2)+minFactor*(maxX(2)-minX(2))/rectFactor 0]' ];
        %vub = [vub ; [minX(1)+maxFactor*(maxX(1)-minX(1))/rectFactor minX(2)+maxFactor*(maxX(2)-minX(2))/rectFactor 0]' ];
        
        vlb = [vlb ; [meanX(1)-meanFactor*(maxX(1)-minX(1)) meanX(2)-meanFactor*(maxX(2)-minX(2)) 0]' ];
        vub = [vub ; [meanX(1)+meanFactor*(maxX(1)-minX(1)) meanX(2)+meanFactor*(maxX(2)-minX(2)) 0]' ];
    end
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

[~,~,Q] = EllipsoidMatricesGenerate(ellParams);
vectorLengths = diag(x'*Q*x);
f = sqrt(mean((vectorLengths-1).^2));
end
