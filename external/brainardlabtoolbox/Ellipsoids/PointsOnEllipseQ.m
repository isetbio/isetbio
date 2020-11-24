function ellPoints = PointsOnEllipseQ(Q,theDirs)
% Generate points on ellipse from Q matrix and point directions
%
% Syntax:
%    ellPoints = PointsOnEllipseQ(Q,theDirs)
%
% Description:
%    Take a set of directions in the plane, and scale so that they lie on the 
%    ellipse defined by matrix Q. These are the set of points x such that 
%        x'*Q*x = 1
%
%    The ellipse is centered at the origin.  
%
%    If Q is passed as a cell array, ellipses for each Q are generated and
%    returned in a cell array.
%
% Inputs:
%    Q               - 2 by 2 positive definite symmetric matrix that defines the
%                      ellipse. This can also be a cell array of such
%                      Qs.
%    theDirs         - 2 by nPoints matrix that define directions in the
%                      plane. These do not have to be normalized.  They
%                      can't be a cell array.
%
% Outputs:
%    ellPoints       - 2 by nPoints matrix of points on the ellipse. If Q is passed
%                      as a cell array, this is a cell array of corresponding ellipse
%                      points.
%
% See also: FitEllipseQ, EllipsoidMatricesGenerate, UnitCircleGenerate

% History:
%   11/23/20  dhb  Got this working, added comments.

% Examples:
%{
    ellParams = [0.5 2 45];
    [~,~,Q] = EllipsoidMatricesGenerate(ellParams,'dimension',2);
    theDirs = UnitCircleGenerate(100);
    ellPoints = PointsOnEllipseQ(Q,theDirs);
    figure; clf; hold on
    plot(ellPoints(1,:),ellPoints(2,:),'ro','MarkerFaceColor','r');
    xlim([-2 2]); ylim([-2 2]); axis('square');
%}


% Normalize each point, compute desired length to reach ellipse, and scale
% accordingly.
if (iscell(Q))
    for cc = 1:length(Q)
        for ii = 1:size(theDirs,2)
            theDirsUse = theDirs(:,ii)/norm(theDirs(:,ii));
            length2 = theDirsUse'*Q{cc}*theDirsUse;
            ellPoints{cc}(:,ii) = theDirsUse/sqrt(length2);
        end
    end
else
    for ii = 1:size(theDirs,2)
        theDirsUse = theDirs(:,ii)/norm(theDirs(:,ii));
        length2 = theDirsUse'*Q*theDirsUse;
        ellPoints(:,ii) = theDirsUse/sqrt(length2);
    end
end

end