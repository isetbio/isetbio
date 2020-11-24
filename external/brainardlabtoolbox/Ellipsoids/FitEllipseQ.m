function [ellParams,AConstraint,Ainv,Q] = FitEllipseQ(theData,varargin)
% Fit an ellipse to points
%
% Syntax:
%    [ellParams,AConstraint,Ainv,Q] = FitEllipseQ(theData)
%
% Description:
%    Fit an ellipse to points.  The fit error minimized is the RMSE
%    of radial length difference in the directions of each data point.
%    This is not a coordinate space invariant method.
%
%    The ellipse is centered at the origin.  
%
%    If data are passed as a cell array, ellipses of same shape are fit to
%    each dataset in the cell array, but their scale is allowed to differ.
%
% Inputs:
%    theData         - Data to be fit, in a 2 by nPoints matrix with each column
%                      one point.  This can also be a cell array of such
%                      points, in which case an ellipse of same shape is
%                      fit to each set.  Return values are also cell arrays
%                      in this case.
%
% Outputs:
%    ellParams       - 3 by one vector of parameters.  First parameter is
%                      reciprocal of major axis length, second is
%                      reciprocol of minor axis length, and third is angle
%                      in degrees of major axis, clockwise from the x-axis.
%                      The lengths of the axis are the equivalent of the
%                      radius of a circle - multiply by 2 to get what would
%                      more intuitively be called the length of the axis.
%                      If more than one ellipse is fit, this and subsequent
%                      outputs come back as a cell array, one for each data
%                      set fit.
%    A               - The ellispe A matrix.  See EllipsoidMatricesGenerate.
%    Ainv            - The inverse of A.
%    Q               - The ellipse Q matrix. See EllipsoidMatricesGenerate. 
%
% Optional key/value pairs:
%    lockAngleAt0    - Lock angle of ellipse at 0.  In this case, the
%                      major/minor axis distinction is relaxed - either x
%                      or y axis of ellipse can be longer.
%    ratioMax        - Maximum major/minor axis ratio.  Default 100.
%    errorScalar     - Multiply RMSE by this in error function, to bring it
%                      into scale where fmincon is happy.  Default 1000.
%                      Try adjusting if search is getting stuck.
%
% See also: PointsOnEllipseQ, EllipsoidMatricesGenerate, UnitCircleGenerate.

% History:
%   11/23/20  dhb  Got this working, added comments.

% Parse input
p = inputParser; p.KeepUnmatched = false;
p.addParameter('lockAngleAt0',false,@islogical);
p.addParameter('ratioMax',100,@isnumeric);
p.addParameter('errorScalar',1000,@isnumeric);
p.parse(varargin{:});

% Set up fmincon
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off');
ratioMax = p.Results.ratioMax;

% If angle locked at 0, either axis can be longer.
if (p.Results.lockAngleAt0)
    angle0 = 0;
    x0 = [1 0.98 angle0];
    vlbVec = [1/ratioMax 1/ratioMax angle0];
    vubVec = [ratioMax ratioMax angle0];
    AConstraint = [];
    bConstraint = [];
    
% If angle is free, enforce major axis longer.  Length (radius) is
% reciprocal of parameter, so we enforce major axis parameter smaller.
else
    angle0 = 0;
    x0 = [0.99 1.01 angle0];
    vlbVec = [1/ratioMax 1/ratioMax -90];
    vubVec = [ratioMax ratioMax 90];
    AConstraint = [1 -1 0];
    bConstraint = 0;
end

% If data comes as a cell array, we add additional parameters that scale
% the ellipse for each data set past the first one.  The convention is
% understood by the fit error routine.
if (iscell(theData))
    x0Cell = x0;
    vlbVecCell = vlbVec;
    vubVecCell = vubVec;
    AConstraintCell = AConstraint;
    for cc = 1:length(theData)
        x0Cell = [x0Cell 1];
        vlbVecCell = [vlbVecCell 1/ratioMax];
        vubVecCell = [vubVecCell ratioMax];
        if (~isempty(AConstraintCell))
            AConstraintCell = [AConstraintCell 0];
        end
    end
    ellParamsRaw = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),x0Cell,AConstraintCell,bConstraint,[],[],vlbVecCell,vubVecCell,[],options);
    
    % Unpack the parameters into cell arrays.
    for cc = 1:length(theData)
        ellParams{cc} = ellParamsRaw(1:3);
        if (cc > 1)
            ellParams{cc}(1:2) = ellParams{cc}(1:2)/ellParamsRaw(3+cc-1);
        end
        
        % Convert parameters to the standard set of ellipse matrices.
        [Atemp, AinvTemp,QTemp] = EllipsoidMatricesGenerate(ellParams{cc},'dimension',2);
        A{cc} = Atemp;
        Ainv{cc} = AinvTemp; 
        Q{cc} = QTemp;
    end
 
% Just one data set.  Fit and convert.  No muss, no fuss.
else
    ellParams = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),x0,AConstraint,bConstraint,[],[],vlbVec,vubVec,[],options);
    [AConstraint,Ainv,Q] = EllipsoidMatricesGenerate(ellParams,'dimension',2);
end

end

% Error function for FitEllipseQ. Internal to this function.
function f = fitErrorFunction(x,theData,errorScalar)

f = 0;
pointCounter = 0;

% If data are a cell array, need to unpack the parameters and compute fit
% to each data set in the cell array.
if (iscell(theData))
    for cc = 1:length(theData)
        xUse = x(1:3);
        if (cc > 1)
            xUse(1:2) = xUse(1:2)/x(3+cc-1);
        end
        [~,~,Q] = EllipsoidMatricesGenerate(xUse,'dimension',2);
        for ii = 1:size(theData{cc},2)
            theDir = theData{cc}(:,ii)/norm(theData{cc}(:,ii));
            length2 = theDir'*Q*theDir;
            ellData = theDir/sqrt(length2);
            f = f + norm(theData{cc}(:,ii)-ellData);
            pointCounter = pointCounter + 1;
        end
    end
    
% If one dataset just do that.
else  
    [~,~,Q] = EllipsoidMatricesGenerate(x,'dimension',2);
    for ii = 1:size(theData,2)
        theDir = theData(:,ii)/norm(theData(:,ii));
        length2 = theDir'*Q*theDir;
        ellData = theDir/sqrt(length2);
        f = f + norm(theData(:,ii)-ellData);
        pointCounter = pointCounter + 1;
    end
    
end

% Convert SSE to RMSE
f = errorScalar*sqrt(f/pointCounter);

end