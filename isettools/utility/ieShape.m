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
%   'circle':
%       center, radius, nSamp
%   'rectangle'
%      rect
%   'line'
%      lineX, lineY
%
%Example:
%
%  vcNewGraphWin; 
%  [h, pts] = ieShape('circle','center',[10 10],'radius', 1,'nSamp',50,'color','r--');
%  axis equal; grid on
%
%  [h, pts] = ieShape('rectangle','rect',[10 10 9 9],'color','b');
%  [h, pts] = ieShape('line','lineX',[0 10],'lineY',[1 3],'color','g');
%
% BW ISETBIO Team

%%
p = inputParser;

p.addRequired('shape',@isstr);
p.addParameter('nSamp',200,@isnumeric);
p.addParameter('center',[0 0],@ismatrix)
p.addParameter('radius',1,@isnumeric);
p.addParameter('rect',[-1 -1 2 2],@isvector)
p.addParameter('lineX',[0 1],@isvector);
p.addParameter('lineY',[0 1],@isvector);
p.addParameter('color','k',@ischar);

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


