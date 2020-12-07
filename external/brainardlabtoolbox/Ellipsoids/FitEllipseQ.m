function [ellParams,AConstraint,Ainv,Q,fitErr] = FitEllipseQ(theData,varargin)
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
%    fitErr          - The minimized fit error.
%
% Optional key/value pairs:
%    lockAngleAt0    - Lock angle of ellipse at 0.
%    ratioMax        - Maximum major/minor axis ratio.  Default 100.
%    errorScalar     - Multiply RMSE by this in error function, to bring it
%                      into scale where fmincon is happy.  Default 1000.
%                      Try adjusting if search is getting stuck because
%                      input data is not close to order unity.
%    initialParams   - Initial parameters. 3 by 1 row vector.  Default
%                      [0.99 1.01 0].
%
% See also: PointsOnEllipseQ, EllipsoidMatricesGenerate, UnitCircleGenerate.

% History:
%   11/23/20  dhb  Got this working, added comments.
%   11/25/20  dhb  Generalized a bit more.

%% Note
%
% At its heart, this routine is simpler than it looks. It's long because it
% brute forces a number of fussy details relating the fact that it can
% accept a cell array of data that can be fit with a set of ellipses of the
% same shape, because there is the option to fit with an ellipse whose
% shape is constrained to align with the axes, and because sometimes the
% full fit works better if initialzed with the result of the aligned fit.
% A bit more work would probably simplify the code, but life seems a little
% too short for that right now.

% Parse input
p = inputParser; p.KeepUnmatched = false;
p.addParameter('lockAngleAt0',false,@islogical);
p.addParameter('ratioMax',100,@isnumeric);
p.addParameter('errorScalar',1000,@isnumeric);
p.addParameter('initialParams',[0.99 1.01 0],@isnumeric);
p.parse(varargin{:});

% Set up fmincon
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off','LargeScale','off');
ratioMax = p.Results.ratioMax;

% If angle locked at 0, either axis can be longer, so no inequality
% constraint.
x0 = p.Results.initialParams;
angle0 = x0(3);
if (p.Results.lockAngleAt0)
    angle0 = 0;
    vlbVec = [1/ratioMax 1/ratioMax angle0];
    vubVec = [ratioMax ratioMax angle0];
    AConstraint = [];
    bConstraint = [];
    
% If angle is free, enforce major axis longer.  Length (radius) is
% reciprocal of parameter, so we enforce major axis parameter smaller.
else
    vlbVec = [1/ratioMax 1/ratioMax -90.2];
    vubVec = [ratioMax ratioMax 90.2];
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
        if (~p.Results.lockAngleAt0)
            AConstraintCell = [AConstraintCell 0];
        end
    end
    
    % If we're not locking angle, do an initial fit with locked angle to
    % help find good initial parameters.
    if (~p.Results.lockAngleAt0)
        x0Temp = x0Cell;
        vlbVecTemp = vlbVecCell;
        vubVecTemp = vubVecCell;
        x0Temp(3) = angle0;
        vlbVecTemp(3) = angle0;
        vubVecTemp(3) = angle0;
        ellParamsRaw0 = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),x0Temp,[],[],[],[],vlbVecTemp,vubVecTemp,[],options);
        ellParamsRaw0 = ellParamsCanonical(ellParamsRaw0);
        fitErr0 = fitErrorFunction(ellParamsRaw0,theData,p.Results.errorScalar);
        
        % Fit with angle locked version as initial parameters.
        %
        % If this fit goes sideways, retain ellParams0 as ellParams2;
        ellParamsRaw2 = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),ellParamsRaw0,AConstraintCell,bConstraint,[],[],vlbVecCell,vubVecCell,[],options);
        ellParamsRaw2 = ellParamsCanonical(ellParamsRaw2);
        fitErr2 = fitErrorFunction(ellParamsRaw2,theData,p.Results.errorScalar);
        if (fitErr0 < fitErr2)
            ellParamsRaw2 = ellParamsRaw0;
            fitErr2 = fitErr0;
        end 
    end
    
    % Fit with initial parameters
    ellParamsRaw1 = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),x0Cell,AConstraintCell,bConstraint,[],[],vlbVecCell,vubVecCell,[],options);
    ellParamsRaw1 = ellParamsCanonical(ellParamsRaw1);
    fitErr1 = fitErrorFunction(ellParamsRaw1,theData,p.Results.errorScalar);

    % If we weren't locked, choose better of passed initial parameters and
    % result of angle locked fit initialization.  Otherwise just happy with
    % what we got from the straight unlocked fit we ran.
    if (~p.Results.lockAngleAt0)
        if (fitErr1 < fitErr2)
            ellParamsRaw = ellParamsRaw1;
            fitErr = fitErr1;
        else
            ellParamsRaw = ellParams2;
            fitErr = fitErr2;
        end
    else
        ellParamsRaw = ellParamsRaw1;
        fitErr = fitErr1;
    end

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
 
% Just one data set.  Fit and convert to canonical form.
else
    % If we're not locking angle, do an initial fit with locked angle to
    % help find good initial parameters.
    if (~p.Results.lockAngleAt0)
        x0Temp = x0;
        vlbVecTemp = vlbVec;
        vubVecTemp = vubVec;
        x0Temp(3) = angle0;
        vlbVecTemp(3) = angle0;
        vubVecTemp(3) = angle0;
        ellParams0 = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),x0Temp,[],[],[],[],vlbVecTemp,vubVecTemp,[],options);
        ellParams0 = ellParamsCanonical(ellParams0);
        fitErr0 = fitErrorFunction(ellParams0,theData,p.Results.errorScalar);
        
        % Fit with angle locked version as initial parameters.
        %
        % If this fit goes sideways, retain ellParams0 as ellParams2;
        ellParams2 = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),ellParams0,AConstraint,bConstraint,[],[],vlbVec,vubVec,[],options);
        ellParams2 = ellParamsCanonical(ellParams2);
        fitErr2 = fitErrorFunction(ellParams2,theData,p.Results.errorScalar);
        if (fitErr0 < fitErr2)
            ellParams2 = ellParams0;
            fitErr2 = fitErr0;
        end 
    end
    
    % Fit with initial parameters
    ellParams1 = fmincon(@(x)fitErrorFunction(x,theData,p.Results.errorScalar),x0,AConstraint,bConstraint,[],[],vlbVec,vubVec,[],options);
    ellParams1 = ellParamsCanonical(ellParams1);
    fitErr1 = fitErrorFunction(ellParams1,theData,p.Results.errorScalar);
    
    % If we weren't locked, choose better of passed initial parameters and
    % result of angle locked fit initialization.  Otherwise just happy with
    % what we got from the straight unlocked fit we ran.
    if (~p.Results.lockAngleAt0)
        if (fitErr1 < fitErr2)
            ellParams = ellParams1;
            fitErr = fitErr1;
        else
            ellParams = ellParams2;
            fitErr = fitErr2;
        end
    else
        ellParams = ellParams1;
        fitErr = fitErr1;
    end
    
    % Convert parameters to the standard set of matrices
    [AConstraint,Ainv,Q] = EllipsoidMatricesGenerate(ellParams,'dimension',2);
end

end

%% Error function for FitEllipseQ. Internal to this function.
function f = fitErrorFunction(x,theData,errorScalar)

% Initialize accumulators
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

%% Convert SSE to RMSE
f = errorScalar*sqrt(f/pointCounter);

end

% Put ellipse parameters into canonical form, with first parameter
% describing the major axis and the angle between -90 and 90
function ellParams = ellParamsCanonical(ellParams)

% Put solution into form with major axis longer than minor axis.  Axis
% length is proportional to reciprocol of parameter.  This is done by
% swapping axis length parameters and rotating by 90 degrees.
if (ellParams(1) > ellParams(2))
    temp = ellParams(1); ellParams(1) = ellParams(2); ellParams(2) = temp;
    ellParams(3) = ellParams(3) + 90;   
end

% Put angle into desired range.  Assumes we were in nearly in range -90 to
% 90 before we added 90 in case of flip.
if (ellParams(3) > 90)
    ellParams(3) = ellParams(3) - 180;
end
if (ellParams(3) < -90)
    ellParams(3) = ellParams(3) + 180;
end

end