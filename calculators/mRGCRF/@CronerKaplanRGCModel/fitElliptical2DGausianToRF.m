function [fittedParams,  rfFunction] = fitElliptical2DGausianToRF(functionName, X, Y, RF, deltaX, minSigma, center)
  
    xyData(:,1) = X;
    xyData(:,2) = Y;
    minX = min(X(:));
    minY = min(Y(:));
    maxX = max(X(:));
    maxY = max(Y(:));
    xRange = maxX-minX;
    yRange = maxY-minY;
    xyRange = max([xRange yRange]);

    assert(ismember(functionName, {'circular Gaussian','elliptical Gaussian'}),...
        sprintf('fitElliptical2DGausianToRF:: unknown function to fit'));
    
    if (strcmp(functionName, 'circular Gaussian'))
        initialParams(1) = 1;  % amplitude
        initialParams(2) = center(1); % x0
        initialParams(3) = center(2); % y0
        initialParams(4) = 0.1; % xSigma
        initialParams(5) = 4; % exponent; 
        fprintf('Fitting a circular Gaussian');
        rfFunction = @circular2DGaussian;
        % Ranges   GAIN        xo                      yo                    sigmaX     exponent
        lb =     [0       center(1)-1000*deltaX     center(2)-1000*deltaX   minSigma       2];
        ub =     [1       center(1)+1000*deltaX     center(2)+1000*deltaX   xyRange       100];
    elseif(strcmp(functionName, 'elliptical Gaussian'))
        initialParams(1) = 1;  % amplitude
        initialParams(2) = center(1); % x0
        initialParams(3) = center(2); % y0
        initialParams(4) = 0.1; % xSigma
        initialParams(5) = 0.1; % ySigma
        initialParams(6) = 20; % rotation (degs)
        initialParams(7) = 4; % exponent; 
        fprintf('Fitting an elliptical Gaussian');
        rfFunction = @elliptical2DGaussian;
        % Ranges   GAIN       xo                 yo                    sigmaX     sigmaY  rotation   exponent
        lb =     [0      center(1)-1000*deltaX   center(2)-1000*deltaX   minSigma   minSigma -360       2.0];
        ub =     [1      center(1)+1000*deltaX   center(2)+1000*deltaX   xyRange    xyRange   360       100.0];
    end
    
    [fittedParams,resnorm,residual,exitflag] = lsqcurvefit(rfFunction, initialParams, xyData, RF,lb,ub);

    doMultiStartSearch = ~true;
    if (doMultiStartSearch)
        problem = createOptimProblem('lsqcurvefit',...
            'x0',fittedParams, ...
            'objective',rfFunction,...
            'lb',lb, ...
            'ub',ub,...
            'xdata',xyData,...
            'ydata',RF);

        displayProgress = 'off'; % 'iter';
        ms = MultiStart(...
            'Display', displayProgress, ...
            'FunctionTolerance', 2e-4, ...
            'UseParallel', true);

        [xmulti,errormulti] = run(ms,problem,20);

        fittedParams = xmulti;
    end
    
end

function f = circular2DGaussian(params,xydata)
    amplitude = params(1);
    xCenter = params(2);
    yCenter = params(3);
    xSigma = params(4);
    exponent = params(5);
    
    xx = xydata(:,1);
    yy = xydata(:,2);
    xx0 = xCenter;
    yy0 = yCenter;

    f = exp( -(abs(xx-xx0)/xSigma).^exponent) .* exp( -(abs(yy-yy0)/xSigma).^exponent);
    f = amplitude * f / max(f(:));
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

    f = exp( -(abs(xx-xx0)/xSigma).^exponent) .* exp( -(abs(yy-yy0)/ySigma).^exponent);
    f = amplitude * f / max(f(:));
    
end