function [DoGparams, theFittedSTF] = genericDifferentOfGaussiansToCompositeSTF(...
    DoGparams, spatialFrequencySupportCPDtoFit, theSTFtoFit, axDoGFitToCompositeSTF, normFactor, multiStartsNum)

     % The DoG model in the frequency domain
     DoGSTF = @(params,spatialFrequency)(...
                    abs(params(1) * ( pi * params(4)^2 * exp(-(pi*params(4)*spatialFrequency).^2) ) - ...
                    params(1)*params(2)/(params(3))^2 * ( pi * (params(4)*params(3))^2 * exp(-(pi*params(4)*params(3)*spatialFrequency).^2) )));
        
     % The optimization objective
     weights = 1 + 0./spatialFrequencySupportCPDtoFit;
     objective = @(p) sum((weights .* (DoGSTF(p, spatialFrequencySupportCPDtoFit) - theSTFtoFit)).^2);
        
     % Multi-start
     problem = createOptimProblem('fmincon',...
          'objective', objective, ...
          'x0', DoGparams.initialValues, ...
          'lb', DoGparams.lowerBounds, ...
          'ub', DoGparams.upperBounds, ...
          'options', optimoptions(...
		 	'fmincon',...
		    'Display', 'off', ...
		    'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
		    'GradObj', 'off', ...
		    'DerivativeCheck', 'off', ...
		    'MaxFunEvals', 10^5, ...
		    'MaxIter', 10^3) ...
          );
      
     ms = MultiStart(...
          'Display', 'off', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);
      
    % Run the multi-start
    DoGparams.finalValues = run(ms, problem, multiStartsNum);
    
    % Visualize the fitted model params and their range
    if (~isempty(axDoGFitToCompositeSTF))
    	     RGCMosaicConstructor.visualize.fittedModelParams(axDoGFitToCompositeSTF, DoGparams, mfilename);
    end


    % Account for normalizationFactor
    DoGparams.finalValues(1) = DoGparams.finalValues(1) * normFactor;

    theFittedSTF.compositeSTF = DoGSTF(DoGparams.finalValues, spatialFrequencySupportCPDtoFit);
    theFittedSTF.centerSTF = DoGparams.finalValues(1) * ( pi * DoGparams.finalValues(4)^2 * exp(-(pi*DoGparams.finalValues(4)*spatialFrequencySupportCPDtoFit).^2) );
    theFittedSTF.surroundSTF = DoGparams.finalValues(1)*DoGparams.finalValues(2)/(DoGparams.finalValues(3))^2  * ( pi * (DoGparams.finalValues(4)*DoGparams.finalValues(3))^2 * exp(-(pi*DoGparams.finalValues(4)*DoGparams.finalValues(3)*spatialFrequencySupportCPDtoFit).^2) );


    differentials = (theSTFtoFit(:) - theFittedSTF.compositeSTF(:));
    theRMSerror = sqrt(sum(differentials.^2)/numel(theSTFtoFit));

    theFittedSTF.RMSE = theRMSerror;
    sfHiRes = logspace(log10(min(spatialFrequencySupportCPDtoFit)), log10(max(spatialFrequencySupportCPDtoFit)), 64);

    theFittedSTF.sfHiRes = sfHiRes;
    theFittedSTF.compositeSTFHiRes = DoGSTF(DoGparams.finalValues, sfHiRes);
    theFittedSTF.centerSTFHiRes = DoGparams.finalValues(1) * ( pi * DoGparams.finalValues(4)^2 * exp(-(pi*DoGparams.finalValues(4)*sfHiRes).^2) );
    theFittedSTF.surroundSTFHiRes = DoGparams.finalValues(1)*DoGparams.finalValues(2)/(DoGparams.finalValues(3))^2 * ( pi * (DoGparams.finalValues(4)*DoGparams.finalValues(3))^2 * exp(-(pi*DoGparams.finalValues(4)*DoGparams.finalValues(3)*sfHiRes).^2) );
end