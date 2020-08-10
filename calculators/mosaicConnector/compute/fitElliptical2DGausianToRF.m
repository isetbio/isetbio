function [fittedParams,  rfFunction] = ...
    fitElliptical2DGausianToRF(X, Y, RF, deltaX, minSigma, inputConeIndices, center)
  
    xyData(:,1) = X;
    xyData(:,2) = Y;
    minX = min(X(:));
    minY = min(Y(:));
    maxX = max(X(:));
    maxY = max(Y(:));
    xRange = maxX-minX;
    yRange = maxY-minY;
    xyRange = max([xRange yRange]);

    if ((~isempty(inputConeIndices)) && (numel(inputConeIndices) == 1))
        initialParams(1) = 1;  % amplitude
        initialParams(2) = center(1); % x0
        initialParams(3) = center(2); % y0
        initialParams(4) = 0.1; % xSigma
        fprintf('Fitting a circular Gaussian');
        rfFunction = @circular2DGaussian;
        % Ranges   GAIN .          xo                      yo               sigmaX    
        lb =     [0.0001    center(1)-100*deltaX   center(2)-100*deltaX   minSigma];
        ub =     [1.0000    center(1)+100*deltaX   center(2)+100*deltaX   xyRange];
    else
        initialParams(1) = 1;  % amplitude
        initialParams(2) = center(1); % x0
        initialParams(3) = center(2); % y0
        initialParams(4) = 0.1; % xSigma
        initialParams(5) = 0.1; % ySigma
        initialParams(6) = 20; % rotation (degs)
        fprintf('Fitting an elliptical Gaussian');
        rfFunction = @elliptical2DGaussian;
        % Ranges   GAIN .          xo                      yo             sigmaX     sigmaY  rotation
        lb =     [0.0001    center(1)-100*deltaX   center(2)-100*deltaX   minSigma   minSigma -360];
        ub =     [1.0000    center(1)+100*deltaX   center(2)+100*deltaX   xyRange    xyRange   360];
    end
    
    [fittedParams,resnorm,residual,exitflag] = lsqcurvefit(rfFunction, initialParams, xyData, RF,lb,ub);
    
end

function f = circular2DGaussian(params,xydata)
    amplitude = params(1);
    xCenter = params(2);
    yCenter = params(3);
    xSigma = params(4);
    
    xx = xydata(:,1);
    yy = xydata(:,2);
    xx0 = xCenter;
    yy0 = yCenter;

    f = amplitude * exp( -((xx-xx0).^2/(xSigma^2) + (yy-yy0).^2/(xSigma^2)) );
end

function f = elliptical2DGaussian(params,xydata)
    amplitude = params(1);
    xCenter = params(2);
    yCenter = params(3);
    xSigma = params(4);
    ySigma = params(5);
    theta = params(6);
    
    xx = xydata(:,1)*cosd(theta) - xydata(:,2)*sind(theta);
    yy = xydata(:,1)*sind(theta) + xydata(:,2)*cosd(theta);
    xx0 = xCenter*cosd(theta) - yCenter*sind(theta);
    yy0 = xCenter*sind(theta) + yCenter*cosd(theta);

    f = amplitude * exp( -((xx-xx0).^2/(xSigma^2) + (yy-yy0).^2/(ySigma^2)) );
end