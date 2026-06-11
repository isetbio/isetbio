function [differenceOfGaussiansModelParams, fittedSpatialTransferFunctions] = spatialTransferFunction(...
        spatialFrequencySupport, spatialTransferFunctions, hiresSFs)
    
    modelFunction = @(params,sf)(...
                    params(1)           * ( pi * (params(2))^2 * exp(-(pi*params(2)*sf).^2) ) - ...
                    params(1)*params(3) * ( pi * (params(2)*params(4))^2 * exp(-(pi*params(2)*params(4)*sf).^2) ) ...
            ); 
        
    %                 Kc   Rc     kS/kC  Rs/Rc
    initialParams = [300   0.05    1e-2   7];
         
    % Upper and lower values of DoG params
    %               Kc       Rc     kS/kC    Rs/Rc  
    lowerBounds   = [1e0     0.0001 1e-4      2  ];
    upperBounds   = [1e10   2.0       1e0       20 ];

    oldoptions = optimoptions('lsqcurvefit');
    options = optimoptions(oldoptions,'MaxFunctionEvaluations',5000, 'FunctionTolerance', 1e-8, 'MaxIterations', 10000);
    
    rgcsNum = size(spatialTransferFunctions,1);
    fittedSpatialTransferFunctions = zeros(rgcsNum, numel(hiresSFs));
    
    % Fit spatial frequency curves for each RGC
    for iRGC = 1:rgcsNum
        
        fprintf('Fitting spatial transfer function %d/%d\n', iRGC, rgcsNum);
        % Fit the model to the data
        fittedParams = lsqcurvefit(modelFunction, initialParams, ...
            spatialFrequencySupport, squeeze(spatialTransferFunctions(iRGC,:)),...
            lowerBounds,upperBounds, options);
       
        % Do a multi-start to find global minimum
        problem = createOptimProblem('lsqcurvefit',...
            'x0',fittedParams, ...
            'objective',modelFunction,...
            'lb',lowerBounds, ...
            'ub',upperBounds,...
            'xdata',spatialFrequencySupport,...
            'ydata',squeeze(spatialTransferFunctions(iRGC,:)));

        displayProgress = 'off'; % 'iter';
        ms = MultiStart(...
            'Display', displayProgress, ...
            'FunctionTolerance', 2e-4, ...
            'UseParallel', true);

        fittedParams = run(ms,problem,50);
        
        % Save fitted params
        differenceOfGaussiansModelParams.centerPeakSensitivities(iRGC) = fittedParams(1);
        differenceOfGaussiansModelParams.centerCharacteristicRadiiDegs(iRGC) = fittedParams(2);
        differenceOfGaussiansModelParams.surroundPeakSensitivities(iRGC) = fittedParams(3)*fittedParams(1);
        differenceOfGaussiansModelParams.surroundCharacteristicRadiiDegs(iRGC) = fittedParams(4)*fittedParams(2);
        
        fprintf('Rc = %f, Ks/Kc = %2.4f, Rs/Rc=%2.1f, S/C int. sens.:%2.2f\n', ...
            fittedParams(2), fittedParams(3), fittedParams(4), fittedParams(3)*fittedParams(4)^2);
        % Compute hires fitted spatial transfer function
        fittedSpatialTransferFunctions(iRGC,:) = modelFunction(fittedParams, hiresSFs);
    end % iRGC
end

