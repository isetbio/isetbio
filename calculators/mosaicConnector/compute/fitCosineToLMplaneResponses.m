function [xFullData, yFullData, xFit, yFit, xIsoResponseContrast, yIsoResponseContrast, ...
    phaseDegs, peakModulation, baserate] = ...
    fitCosineToLMplaneResponses(LMresponseMagnitude, LMangles, isoResponseLevelPercentOfMax)

    x = LMresponseMagnitude .* cosd(LMangles(:));
    y = LMresponseMagnitude .* sind(LMangles(:));

    % Make data symmetric
    xFullData = cat(2,x', -x');
    yFullData = cat(2,y', -y');
    thetaDegs = cat(2,LMangles, LMangles+180);
    response = cat(2, LMresponseMagnitude', LMresponseMagnitude');
    
    % Fit abs(cos) model
    modelFunction = @(params, thetaDegs)(params(1) + params(2)*abs(cosd(thetaDegs-params(3))));

    % Initial params and ranges
    initialParams = [0 10 45];
    lowerBounds = [0 0 -360];
    upperBounds = [max(response)  max(response)*2 360];
    
    oldoptions = optimoptions('lsqcurvefit');
    options = optimoptions(oldoptions,'MaxFunctionEvaluations',5000, 'FunctionTolerance', 1e-8, 'MaxIterations', 10000);
    fittedParams = lsqcurvefit(modelFunction, initialParams, thetaDegs, response, lowerBounds,upperBounds, options);
    
    % Do a multi-start to find global minimum
    problem = createOptimProblem('lsqcurvefit',...
            'x0',fittedParams, ...
            'objective',modelFunction,...
            'lb',lowerBounds, ...
            'ub',upperBounds,...
            'xdata',thetaDegs,...
            'ydata',response);

    displayProgress = 'off'; % 'iter';
    ms = MultiStart(...
            'Display', displayProgress, ...
            'FunctionTolerance', 2e-4, ...
            'UseParallel', true);

   [fittedParams,errormulti] = run(ms,problem,50);
        
    
   baserate = fittedParams(1);
   peakModulation = fittedParams(2);
   phaseDegs = fittedParams(3);
    
   angleFit = -90:2:268;
   responseFit = modelFunction(fittedParams, angleFit);
 
   xFit = (responseFit .* cosd(angleFit));
   yFit = (responseFit .* sind(angleFit));    
    
   % Find iso-response contrast for some response level
   isoResponseLevel = isoResponseLevelPercentOfMax*max(responseFit);
   constrastBoostToReachIsoResponseLevel = (isoResponseLevel-baserate)./(responseFit-baserate);
    
   gain = max(responseFit);
   xIsoResponseContrast = gain * constrastBoostToReachIsoResponseLevel .* cosd(angleFit);
   yIsoResponseContrast = gain * constrastBoostToReachIsoResponseLevel .* sind(angleFit);
end

