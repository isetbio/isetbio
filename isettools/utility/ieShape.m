function [h, pts] = ieShape(shape, varargin)
% Draw a shape on the current window
%
% Syntax:
%   [h, pts] = ieShape(type, [varargin])
%
% Description:
%    This function will draw a shape on the current window.
%
%    Examples in the code.
%
% Inputs:
%    shape     - Required input variable, with possible options of 'line', 
%                'circle', 'rectangle', or 'ellipse' The shapes have the
%                following optional parameters listed below.
%
% Outputs:
%    h         - Handle containin the shape's information
%    pts       - collection of points to draw
%
% Optional key/value pairs:
%    'color'   - Optional colorspec (e.g. 'b', 'g') information for shapes.
%                Default is black (k).
%    'center'  - Optional key for ellipses and circles. Default is [0 0].
%    'ellipseParameters'
%              - Optional key for ellipses, that takes the values major, 
%                minor, and angle in an array. Default is [1 1 0].
%    'radius'  - Optional key for circles, contains radius. Default is 1.
%    'nSamp'   - Optional key for circles, contains the number of sample
%                points that are used to make up the shape. Default is 20.
%    'rect'    - Optional key for rectangles, contains starting x and y, 
%                and then the length for x and y. Default is  [-1 -1 2 2].
%    'lineX'   - Optional key for lines for the starting X coordinate and
%                the length of the line along the x-axis. Default is [0 1].
%    'lineY'   - Optional key for lines for starting Y coordinate and the
%                length of the line along the y-axis. Default is [0 1].
%
% Notes:
%    * [Note: BW.  We could add a fill color. Set 'MarkerFaceColor']

% History:
%    xx/xx/xx  BW   ISETBIO Team
%    11/22/17  jnm  Formatting
%    01/19/18  jnm  Formatting update to match wiki, add defaults for OKVP.

% Examples:
%{
    vcNewGraphWin; 
    [h, pts] = ieShape('circle', 'center', [10 10], 'radius', 1, ...
        'nSamp', 50, 'color', 'r--');
    axis equal;
    grid on

    [h, pts] = ieShape('rectangle', 'rect', [10 10 9 9], 'color', 'b');
    [h, pts] = ieShape('line', 'lineX', [0 10], 'lineY', [1 3], ...
        'color', 'g');

    [h, pts] = ieShape('ellipse');
    [h, pts] = ieShape('ellipse', 'ellipseParameters', [1 2 45]);
    [h, pts] = ieShape('ellipse', 'center', [0 0 ; 3 3], ...
        'ellipseParameters', [1 2 32; 2 1 70], 'color', 'b');
%}

%%
p = inputParser;

p.addRequired('shape', @isstr);
p.addParameter('nSamp', 20, @isnumeric);
p.addParameter('center', [0 0], @ismatrix);
p.addParameter('radius', 1, @isnumeric);
p.addParameter('rect', [-1 -1 2 2], @isvector)
p.addParameter('lineX', [0 1], @isvector);
p.addParameter('lineY', [0 1], @isvector);
p.addParameter('color', 'k', @ischar);
p.addParameter('fillColor', 'k', @ischar);
p.addParameter('fillArray', [], @isnumeric);
p.addParameter('ellipseParameters', [1 1 0], @ismatrix);  % y, x, rho

p.parse(shape, varargin{:});

nSamp = p.Results.nSamp;

%%  Switch on shape

switch shape
    case 'circle'
        hold on;
        center = p.Results.center;
        radius = p.Results.radius;
        nCircles = size(center, 1);
        
        % It is OK to send in a single value for the radius
        if isscalar(radius) && nCircles > 1
            radius = repmat(radius, nCircles, 1); 
        end
        % It is OK to send in a single color
        if length(p.Results.color) == 1 && nCircles > 1
            colors = repmat(p.Results.color, nCircles, 1); 
        else
            colors = p.Results.color;
        end
        
        % We want a fill color, too, don't we. Wonder how to do that?
        % Also for multiple circles, keep hold on and make axis equal, of
        % course. Otherwise it's not a circle.
        hold on
        for ii = 1:nCircles
            pts = circle(center(ii, :), radius(ii), nSamp);
            h = plot(pts(:, 2), pts(:, 1), colors(ii));
        end
        axis equal
        hold off
        
    case 'ellipse'
        hold on;
        center = p.Results.center;
        
        radius = p.Results.radius;
        nEllipses = size(center, 1);
        
        % It is OK to send in a single ellipse parameter. The ellipse will
        % be replicated at all the center locations.
        ellipseParameters = p.Results.ellipseParameters;
        if size(ellipseParameters, 1) == 1 && nEllipses > 1
            ellipseParameters = repmat(ellipseParameters, nEllipses, 1);
        end
        % It is OK to send in a single color
        if length(p.Results.color) == 1 && nEllipses > 1
            colors = repmat(p.Results.color, nEllipses, 1); 
        else
            colors = p.Results.color;
        end
        
        if length(p.Results.fillArray) == 1 && nEllipses > 1 && ...
                ~isempty(p.Results.fillColor)
            fillArray = repmat(p.Results.fillArray, nEllipses, 1);
            
            fillArrayMax = 1;
        else
            fillArray = p.Results.fillArray;
            fillArrayMax = max(fillArray(:));
        end
        % We want a fill color, too, don't we. Wonder how to do that?
        hold on
        ptsCircle = circle([0, 0], radius, nSamp);
        % [xc, yc] = ind2sub(szEllMatrix, ii);
        % szEllMatrix = [sqrt(nEllipses) sqrt(nEllipses)];
        
        
        for ii = 1:nEllipses
            % ePts = pts * eMatrix;
            % ePts = pts * diag(major / minor axis) * ...
            %     rotmat(third parameter)
            % ellipseParameters(ii, 1:2) = sqrt(2) ...
            %     * ellipseParameters(ii, 1:2) ...
            %     ./ norm([ellipseParameters(ii, 1:2)]);
            D = diag(ellipseParameters(ii, 1:2));
            R = [cosd(ellipseParameters(ii, 3)) ...
                -sind(ellipseParameters(ii, 3));
                sind(ellipseParameters(ii, 3)) ...
                cosd(ellipseParameters(ii, 3))];
            pts = bsxfun(@plus, ptsCircle * D * R, center(ii, :));
            %pts = ptsCircle * D * R + repmat(center, nSamp, 1);
            if isempty(fillArray)
                h = plot(pts(:, 2), pts(:, 1), colors(ii), ...
                    'linewidth', .2);
            else
                h = patch(pts(:, 2), pts(:, 1), fillArray(ii) ./ ...
                    fillArrayMax);%, 'linewidth', .2);
            end
        end
        % pts = (radius(ii) * (ellipseMatrix{xc, yc} ./ ...
        %      norm(ellipseMatrix{xc, yc}(:))) * (ptsCircle - ones(200, ...
        %      1) * center(ii, :))')';
        % %([0 1; 1 0]) * (ptsCircle - ones(200, 1) * center(ii, :))';

        axis equal
        hold off

    case 'rectangle'
        % rect = [10 10 50 50];
        % ieDrawShape('rectangle', rect);
        h = rectangle('Position', p.Results.rect);
        set(h, 'edgecolor', p.Results.color);
        pts = p.Results.rect;
    case 'line'
        % X = [0 96];
        % Y = [32 32];
        % ieDrawShape(obj, 'line', X , Y);
        h = line(p.Results.lineX, p.Results.lineY, 'LineWidth', 2, ...
            'color', p.Results.color);
        pts(:, 1) = p.Results.lineX(:);
        pts(:, 2) = p.Results.lineY(:);
    otherwise 
        error('Unknown shape %s\n', shape);
end

end

function pts = circle(center, r, nPts)
% Calculate the circle
%
% Syntax:
%   circle(center, r, nPts);
%
% Description:
%    Using the provided center point, radius, and the number of points to
%    bound the circle, calculate the points.
%
% Inputs:
%    center - Center coordinates for the circle
%    r      - Radius of the circle
%    nPts   - Number of points to bound the circle with
%
% Outputs:
%    pts    - Collection of coordinates for the bounding points
%
% Optional key/value pairs:
%    None.
%
t = linspace(0, 2 * pi, nPts)'; 
pts = zeros(length(t), 2);
pts(:, 1) = r .* cos(t) + center(1); 
pts(:, 2) = r .* sin(t) + center(2); 

end
