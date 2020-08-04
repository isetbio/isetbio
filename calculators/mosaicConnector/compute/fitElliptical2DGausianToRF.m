function [rfPosition, rfSigmas, rfGain, rfRotationDegs, rfExponent, rfFunction] = ...
    fitElliptical2DGausianToRF(X, Y, RF, deltaX, center)

    initialParams(1) = 1;  % amplitude
    initialParams(2) = center(1); % x0
    initialParams(3) = center(2); % y0
    initialParams(4) = 0.1; % xSigma
    initialParams(5) = 0.1; % ySigma
    initialParams(6) = 20; % rotation (degs)
    initialParams(7) = 0.5; % exponent
    xyData(:,1) = X;
    xyData(:,2) = Y;
    minX = min(X(:));
    minY = min(Y(:));
    maxX = max(X(:));
    maxY = max(Y(:));
    xRange = maxX-minX;
    yRange = maxY-minY;
    xyRange = max([xRange yRange]);
    lb = [0.0001 center(1)-100*deltaX center(2)-100*deltaX deltaX   deltaX -360 1.0];
    ub = [1   center(1)+100*deltaX center(2)+100*deltaX xyRange xyRange  360 1.0];
    
    [fittedParams,resnorm,residual,exitflag] = lsqcurvefit(@elliptical2DGaussian, initialParams, xyData, RF,lb,ub);

    rfGain = fittedParams(1);
    rfPosition = [fittedParams(2) fittedParams(3)];
    rfSigmas = [fittedParams(4) fittedParams(5)];
    rfRotationDegs = fittedParams(6);
    rfExponent = fittedParams(7);
    rfFunction = @elliptical2DGaussian;
end

function f = elliptical2DGaussian(params,xydata)
    amplitude = params(1);
    xCenter = params(2);
    yCenter = params(3);
    xSigma = params(4);
    ySigma = params(5);
    theta = params(6);
    exponent = params(7);
    
    xx = xydata(:,1)*cosd(theta) - xydata(:,2)*sind(theta);
    yy = xydata(:,1)*sind(theta) + xydata(:,2)*cosd(theta);
    xx0 = xCenter*cosd(theta) - yCenter*sind(theta);
    yy0 = xCenter*sind(theta) + yCenter*cosd(theta);

    f = exp( -((xx-xx0).^2/(xSigma^2) + (yy-yy0).^2/(ySigma^2)) );
    
    % Flat top
    f = amplitude * f .^ exponent;
    
end