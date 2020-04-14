function  [xEllipse,  yEllipse, semiAxes, rfCenter, noFit] = fitEllipseToContour(xContour,  yContour)
    
    xx = [xContour(:) yContour(:)];
    %[center, a, b, rotationTheta, noFit] = fitellipse(xx, 'linear');  % linear, bookstein contraint
     [center, a, b, rotationTheta, noFit] = fitellipse(xx, 'linear', 'constraint', 'trace');  % Trace contrainst
   %[center, a, b, rotationTheta, noFit] =  fitellipse(xx);  % nonlinear fit
    
    if (noFit)
        semiAxes = [nan nan];
        rfCenter = [nan nan];
        % Return the input points
        xEllipse = xContour;
        yEllipse = yContour;
    else
        % Return the fitted ellipse points
        semiAxes = [a b];
        rfCenter = center;
        [xEllipse,  yEllipse] = makeellipse(center, a, b, rotationTheta);
    end
end

function [x,y] = makeellipse(center, a, b, rotationTheta)
    % form the parameter vector
    npts = 100;
    t = linspace(0, 2*pi, npts);

    % Rotation matrix
    Q = [cos(rotationTheta), -sin(rotationTheta); sin(rotationTheta) cos(rotationTheta)];
    
    % Ellipse points
    X = Q * [a * cos(t); b * sin(t)] + repmat(center, 1, npts);
    
    x = X(1,:);
    y = X(2,:);
    
end
