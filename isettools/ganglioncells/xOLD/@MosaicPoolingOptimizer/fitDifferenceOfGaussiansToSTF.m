function [DoGparams, theFittedSTF] = fitDifferenceOfGaussiansToSTF(...
    spatialFrequencySupportCPD, theSTF, RcDegsInitialEstimate, rangeForRc, multiStartsNum)

    normFactor = max(theSTF(:));
    theSTF = theSTF / normFactor;

    % Only fit the part of the STF that contains the main peak
    [spatialFrequencySupportCPDtoFit, theSTFtoFit] = MosaicPoolingOptimizer.stfPortionToAnalyze(spatialFrequencySupportCPD, theSTF);

    % DoG param initial values and limits: center gain, kc
    Kc = struct(...    
        'low', 0.1, ...
        'high', 1e5, ...
        'initial', 1);

    % DoG param initial values and limits: Ks/Kc ratio
    KsToKc = struct(...
        'low', 1e-4, ...
        'high', 0.2, ...
        'initial', 0.1);

    % DoG param initial values and limits: RsToRc ratio
    RsToRc = struct(...
        'low', 1.5, ...
        'high', 20, ...
        'initial', 5);

    % DoG param initial values and limits: RcDegs
    if (~isempty(rangeForRc))
            RcDegs = struct(...
                'low', rangeForRc(1), ...
                'high', rangeForRc(3), ...
            'initial', rangeForRc(2));
    else
            RcDegs = struct(...
                'low', RcDegsInitialEstimate/sqrt(2.0), ...
                'high', RcDegsInitialEstimate*200, ...
            'initial', RcDegsInitialEstimate*5);
    end

     %                          Kc           kS/kC             RsToRc            RcDegs    
     DoGparams.initialValues = [Kc.initial   KsToKc.initial    RsToRc.initial    RcDegs.initial];
     DoGparams.lowerBounds   = [Kc.low       KsToKc.low        RsToRc.low        RcDegs.low];
     DoGparams.upperBounds   = [Kc.high      KsToKc.high       RsToRc.high       RcDegs.high];
     DoGparams.names         = {'Kc',        'kS/kC',         'RsToRc',         'RcDegs'};
     DoGparams.scaling       = {'log',       'log',           'linear',         'linear'};
     
     % The DoG model in the frequency domain
     DoGSTF = @(params,spatialFrequency)(...
                    abs(params(1)       * ( pi * params(4)^2             * exp(-(pi*params(4)*spatialFrequency).^2) ) - ...
                    params(1)*params(2) * ( pi * (params(4)*params(3))^2 * exp(-(pi*params(4)*params(3)*spatialFrequency).^2) )));
        
     % The optimization objective
     weights = 1 + 0./spatialFrequencySupportCPDtoFit;
     objective = @(p) sum((weights .* (DoGSTF(p, spatialFrequencySupportCPDtoFit) - theSTFtoFit)).^2);

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
          'x0', DoGparams.initialValues, ...
          'lb', DoGparams.lowerBounds, ...
          'ub', DoGparams.upperBounds, ...
          'options', options...
          );
      
     ms = MultiStart(...
          'Display', 'off', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);
      
     % Run the multi-start
     DoGparams.finalValues = run(ms, problem, multiStartsNum);

    
     % Account for normalizationFactor
     DoGparams.finalValues(1) = ...
         DoGparams.finalValues(1) * normFactor;

     theFittedSTF.compositeSTF = DoGSTF(DoGparams.finalValues, spatialFrequencySupportCPD);
     theFittedSTF.centerSTF = DoGparams.finalValues(1) * ( pi * DoGparams.finalValues(4)^2 * exp(-(pi*DoGparams.finalValues(4)*spatialFrequencySupportCPD).^2) );
     theFittedSTF.surroundSTF = DoGparams.finalValues(1)*DoGparams.finalValues(2) * ( pi * (DoGparams.finalValues(4)*DoGparams.finalValues(3))^2 * exp(-(pi*DoGparams.finalValues(4)*DoGparams.finalValues(3)*spatialFrequencySupportCPD).^2) );
     
     sfHiRes = logspace(log10(0.1), log10(100), 64);
     theFittedSTF.sfHiRes = sfHiRes;
     theFittedSTF.compositeSTFHiRes = DoGSTF(DoGparams.finalValues, sfHiRes);
     theFittedSTF.centerSTFHiRes = DoGparams.finalValues(1) * ( pi * DoGparams.finalValues(4)^2 * exp(-(pi*DoGparams.finalValues(4)*sfHiRes).^2) );
     theFittedSTF.surroundSTFHiRes = DoGparams.finalValues(1)*DoGparams.finalValues(2) * ( pi * (DoGparams.finalValues(4)*DoGparams.finalValues(3))^2 * exp(-(pi*DoGparams.finalValues(4)*DoGparams.finalValues(3)*sfHiRes).^2) );

     % figure(1); clf;
     % plot(spatialFrequencySupportCPD, theSTF, 'ks');
     % hold on
     % plot(spatialFrequencySupportCPDtoFit, theSTFtoFit, 'b.');
     % set(gca, 'XScale', 'log', 'XLim', [0.1 100]);
     % plot(theFittedSTF.sfHiRes, theFittedSTF.compositeSTFHiRes, 'r-')   
     % drawnow

     % DoGparams.lowerBounds(1)
     % DoGparams.finalValues(1)
     % DoGparams.upperBounds(1)
     % DoGparams.lowerBounds(2:end)
     % DoGparams.finalValues(2:end)
     % DoGparams.upperBounds(2:end)
     % pause