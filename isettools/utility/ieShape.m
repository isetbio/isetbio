function [h,pts] = ieShape(shape,varargin)
%IESHAPE - Draw a shape on the current window
%
%    [h,pts] = ieShape(type,varargin)
%
% Required input:
%   shape: 'circle', 'rectangle', or 'line'
%
% Optional parameter for all shapes
%   'color'  - Colorspec (e.g., 'b', 'g')
%
% Optional Parameter/Vals for the specific shapes
%   'ellipse'
%      center, ellipseParameters (A,B,rho) (major, minor, angle)
%   'circle':
%       center, radius, nSamp
%   'rectangle'
%      rect
%   'line'
%      lineX, lineY
%
%Examples:
%
%  vcNewGraphWin; 
%  [h, pts] = ieShape('circle','center',[10 10],'radius', 1,'nSamp',50,'color','r--');
%  axis equal; grid on
%
%  [h, pts] = ieShape('rectangle','rect',[10 10 9 9],'color','b');
%  [h, pts] = ieShape('line','lineX',[0 10],'lineY',[1 3],'color','g');
%
%  [h, pts] = ieShape('ellipse');
%  [h, pts] = ieShape('ellipse','ellipseParameters',[1 2 45]);
%  [h, pts] = ieShape('ellipse','center',[0 0 ; 3 3],'ellipseParameters',[1 2 32; 2 1 70],'color','b');
%
% BW ISETBIO Team

%%
p = inputParser;

p.addRequired('shape',@isstr);
p.addParameter('nSamp',20,@isnumeric);
p.addParameter('center',[0 0],@ismatrix);
p.addParameter('radius',1,@isnumeric);
p.addParameter('rect',[-1 -1 2 2],@isvector)
p.addParameter('lineX',[0 1],@isvector);
p.addParameter('lineY',[0 1],@isvector);
p.addParameter('color','k',@ischar);
p.addParameter('ellipseParameters',[1 1 0],@ismatrix);  % y,x,rho

p.parse(shape,varargin{:});

nSamp = p.Results.nSamp;

%%  Switch on shape

switch shape
    case 'circle'
        hold on;
        center = p.Results.center;
        radius = p.Results.radius;
        nCircles = size(center,1);
        
        % It is OK to send in a single value for the radius
        if isscalar(radius) && nCircles > 1
            radius = repmat(radius,nCircles,1); 
        end
        % It is OK to send in a single color
        if length(p.Results.color) == 1 && nCircles > 1
            colors = repmat(p.Results.color,nCircles,1); 
        else
            colors = p.Results.color;
        end
        
        % We want a fill color, too, don't we.  Wonder how to do that?
        % Also for multiple circles, keep hold on and make axis equal, of
        % course.  Otherwise it's not a circle.
        hold on
        for ii=1:nCircles
            pts = circle(center(ii,:),radius(ii),nSamp);
            h = plot(pts(:,2),pts(:,1),colors(ii));
        end
        axis equal
        hold off
        
    case 'ellipse'
        hold on;
        center = p.Results.center;
        
        radius = p.Results.radius;
        nEllipses = size(center,1);
        
        % It is OK to send in a single ellipse parameter.  The ellipse will
        % be replicated at all the center locations.
        ellipseParameters = p.Results.ellipseParameters;
        if size(ellipseParameters,1)==1 && nEllipses > 1
            ellipseParameters = repmat(ellipseParameters,nEllipses,1);
        end
        % It is OK to send in a single color
        if length(p.Results.color) == 1 && nEllipses > 1
            colors = repmat(p.Results.color,nEllipses,1); 
        else
            colors = p.Results.color;
        end
        
        % We want a fill color, too, don't we.  Wonder how to do that?
        % Also for multiple circles, keep hold on and make axis equal, of
        % course.  Otherwise it's not a circle.
        hold on
        ptsCircle = circle([0,0],radius,nSamp);
        % [xc,yc] = ind2sub(szEllMatrix,ii);
        % szEllMatrix = [sqrt(nEllipses) sqrt(nEllipses)];
        for ii=1:nEllipses
            % ePts = pts*eMatrix;
            % ePts = pts * diag(major/minor axis)* rotmat(third parameter)
            D = diag(ellipseParameters(ii,1:2));
            R = [cosd(ellipseParameters(ii,3)) -sind(ellipseParameters(ii,3));
                sind(ellipseParameters(ii,3))   cosd(ellipseParameters(ii,3))];
            pts = bsxfun(@plus,ptsCircle*D*R,center(ii,:));
            %pts = ptsCircle*D*R + repmat(center,nSamp,1);
            h = plot(pts(:,2),pts(:,1),colors(ii),'linewidth',.2);
        end
        % pts = (radius(ii)*(ellipseMatrix{xc,yc}./norm(ellipseMatrix{xc,yc}(:)))*(ptsCircle-ones(200,1)*center(ii,:))')';%([0 1; 1 0])*(ptsCircle-ones(200,1)*center(ii,:))';

        axis equal
        hold off

    case 'rectangle'
        % rect = [10 10 50 50];
        % ieDrawShape('rectangle',rect);
        h = rectangle('Position',p.Results.rect);
        set(h,'edgecolor',p.Results.color);
        pts = p.Results.rect;
    case 'line'
        % X = [0 96]; Y = [32 32];
        % ieDrawShape(obj,'line',X ,Y);
        h = line(p.Results.lineX,p.Results.lineY,'LineWidth',2,'color',p.Results.color);
        pts(:,1) = p.Results.lineX(:);
        pts(:,2) = p.Results.lineY(:);
    otherwise 
        error('Unknown shape %s\n',shape);
end

end

%
function pts = circle(center,r,nPts)

% Calculate the circle
t = linspace(0,2*pi,nPts)'; 
pts = zeros(length(t),2);
pts(:,1) = r.*cos(t) + center(1); 
pts(:,2) = r.*sin(t) + center(2); 

end


