function [h,pts] = ieShape(shape,varargin)
% Draw a shape on the current window
%
%    [h,pts] = ieShape(type,varargin)
%
% Inputs:
%   shape: circle, rectangle, line - other stuff later.
%
% Parameter/Vals for the specific shapes
%   nSamp
%   center, radius  (circle)
%   rect (rectangle)
%   lineX, lineY  (line)
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
% BW ISETBIO Team 2009

%%
p = inputParser;

p.addRequired('shape',@isstr);
p.addParameter('nSamp',200,@isnumeric);
p.addParameter('center',[0 0],@isvector)
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
        % ieDrawShape('circle',[20 20],10);
        hold on;
        center = p.Results.center;
        radius = p.Results.radius;
        pts = circle(center,radius,nSamp);
        
        h = plot(pts(:,2),pts(:,1),p.Results.color);

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


