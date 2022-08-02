function [theRetinalRFparamsStruct, theFittedRF] = fitDiffOfGaussianCenterAndDoubleExponentSurroundModelToRF(spatialSupportDegs, RF, RcDegs)

    %                             Kc   RcDegs    KsToKcPeakRatio  narrowToWideVolumeRatio  RwideDegs  RnarrowToRwideRatio
    DoGDERFparams.initialValues = [  1    0.05    0.3              0.3                      0.2         0.1];
    DoGDERFparams.lowerBounds   = [ 1e-3  0.1/60  0.00             0.01                     0.1/60      0.01];
    DoGDERFparams.upperBounds   = [ 1e2   2       1                1.0                      40          10.0];

    if (~isempty(RcDegs))
        DoGDERFparams.initialValues(2) = RcDegs;
        DoGDERFparams.lowerBounds(2) = RcDegs;
        DoGDERFparams.upperBounds(2) = RcDegs;
    end

    % The optimization objective
    objective = @(p) sum((RetinaToVisualFieldTransformer.diffOfGaussianCenterAndDoubleExponentSurround(p, spatialSupportDegs) - RF).^2, 'all');

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
          'x0', DoGDERFparams.initialValues, ...
          'lb', DoGDERFparams.lowerBounds, ...
          'ub', DoGDERFparams.upperBounds, ...
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
     theFittedRF = RetinaToVisualFieldTransformer.diffOfGaussianCenterAndDoubleExponentSurround(theFittedParamsVector, spatialSupportDegs);

     % Return the fitted params struct
     theRetinalRFparamsStruct = RetinaToVisualFieldTransformer.paramsVectorToParamsStruct(theFittedParamsVector, 'GaussianCenterDoubleExponentSurroundBased');
end