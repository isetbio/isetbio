function [theRetinalRFparamsStruct, theFittedRF] = fitDoGModelToRF(spatialSupportDegs, RF, rotationDegs, RcDegs)  
    %                             Kc   RcDegs   Rs/Rc integratedS/C
    DoGRFparams.initialValues = [  1    0.05    5     0.5];
    DoGRFparams.lowerBounds   = [ 1e-3  0.1/60  1     0.0];
    DoGRFparams.upperBounds   = [ 1e2   2       10     10];

    if (~isempty(RcDegs))
        DoGRFparams.initialValues(2) = RcDegs;
        DoGRFparams.lowerBounds(2) = RcDegs;
        DoGRFparams.upperBounds(2) = RcDegs;
    end

    unrotatedRF = imrotate(RF, -(rotationDegs+90), "bilinear", "crop");

    % The optimization objective
    objective = @(p) sum((RetinaToVisualFieldTransformer.diffOfGaussiansRF(p, spatialSupportDegs) - unrotatedRF).^2, 'all');

    % Ready to fit
    options = optimset(...
            'Display', 'off', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^3);
        
     % Multi-start
     problem = createOptimProblem('fmincon',...
          'objective', objective, ...
          'x0', DoGRFparams.initialValues, ...
          'lb', DoGRFparams.lowerBounds, ...
          'ub', DoGRFparams.upperBounds, ...
          'options', options...
          );
      
     ms = MultiStart(...
          'Display', 'off', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);
      
     % Run the multi-start
     multiStartsNum = 64;
     theFittedParamsVector = run(ms, problem, multiStartsNum);

     % Generate the fitted RF
     theFittedRF = RetinaToVisualFieldTransformer.diffOfGaussiansRF(theFittedParamsVector, spatialSupportDegs);

     % Return the fitted params struct
     theRetinalRFparamsStruct = RetinaToVisualFieldTransformer.paramsVectorToParamsStruct(theFittedParamsVector, 'GaussianCenterGaussianSurroundBased');
end