%
% RGCMosaicConstructor.temporalFilterEngine.impulseResponseFunctionFromTTF
%

function theInnerRetinaTTF = delayHighPassFilterToTransformPhotocurrentTTFtoTargetTTF(...
                temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF, theInnerRetinaTTFfromDirectDeconvolution)

    theTargetTTF = theTargetTTF / max(abs(theTargetTTF));

    useParallel = true;
    multiStartsNum = 32;

    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1100 1000]);
    ax = subplot(1,1,1);

    [~, modelParams.initialValues, ...
        modelParams.lowerBounds, ...
        modelParams.upperBounds, ...
        modelParams.names] = delayLeadLagFilter([],temporalFrequencySupportHz);

    modelParams.finalValues = modelParams.initialValues;
    modelParams.scaling = 'linear';





    hFig = figure(2000); clf;
    set(hFig, 'Position', [500 300 1100 1000]);
    ax = subplot(1,1,1);
    figNo = 2001;

    objectiveFunctionToMinimize = @(x)theObjectiveFunctionToMinimize(x, temporalFrequencySupportHz, ...
        theTargetTTF, thePhotocurrentsBasedTTF, ax, figNo, modelParams);

    
    problem = createOptimProblem('fmincon',...
          'objective', objectiveFunctionToMinimize , ...
          'x0', modelParams.initialValues, ...
          'lb', modelParams.lowerBounds, ...
          'ub', modelParams.upperBounds, ...
          'options', optimoptions(...
            'fmincon',...
            'Display', 'none', ...
            'Algorithm', 'sqp',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^4) ...
          );


    solverType = 'global-search';

     switch (solverType)
            case 'multi-start'
                % Setup the multi-start solver
                ms = MultiStart(...
                'Display', 'iter', ...
                'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
                'UseParallel', useParallel);
    
                % Run the multi-start
                modelParams.finalValues = run(ms, problem, multiStartsNum);
    
            case 'global-search'
    
                gs = GlobalSearch;
                gs = GlobalSearch(gs,'XTolerance',1e-3,'StartPointsToRun','bounds');
                modelParams.finalValues = run(gs,problem);
    
    
            case 'fmincon'
                modelParams.finalValues = fmincon(problem);
    
            otherwise
                error('Uknown solver type: ''%s''.', solverType)
                
     end


    modelParams.finalValues
    theInnerRetinaTTF = delayLeadLagFilter(modelParams.finalValues,temporalFrequencySupportHz);
end


function theResidual = theObjectiveFunctionToMinimize(theCurrentParams, ...
    temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF,  ax, figNo, modelVariables)

    % Composite (delay + leadLagFiter)
    [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = delayLeadLagFilter(theCurrentParams, temporalFrequencySupportHz);

    modelVariables.finalValues = theCurrentParams;

    RGCMosaicConstructor.visualize.fittedModelParams(ax, modelVariables, 'TTF fit');

    
    performFFTshift = true;
    zeroPaddingLength = 512;
    theCurrentInnerRetinaFilterResponseData = RGCMosaicConstructor.temporalFilterEngine.impulseResponseFunctionFromTTF(...
        theCurrentInnerRetinaFilterTTF, temporalFrequencySupportHz, performFFTshift, zeroPaddingLength);

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1000 1000]);
    ax = subplot(2,1,1);
    plot(theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds*1e3, theCurrentInnerRetinaFilterResponseData.amplitude, 'k-');
    set(gca, 'XLim', [500 1500])
    ax = subplot(2,1,2);
    plot(temporalFrequencySupportHz, abs(theTargetTTF), 'r-');
    hold on
    plot(temporalFrequencySupportHz, abs(thePhotocurrentsBasedTTF), 'k-');
    plot(temporalFrequencySupportHz, abs(theCurrentInnerRetinaFilterTTF), 'b-');
    plot(temporalFrequencySupportHz, abs(theCurrentInnerRetinaFilterTTF .* thePhotocurrentsBasedTTF), 'm--');
    set(gca, 'XLim', [0.5 200], 'XScale', 'log', 'XTick', [1 3 10 30 100])

    theResidual = norm(theTargetTTF - theCurrentInnerRetinaFilterTTF.*thePhotocurrentsBasedTTF);
end


function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = delayLeadLagFilter(theCurrentParams, temporalFrequencySupportHz)
    
    % gain
    initialValues(1) = 0.13;
    lowerBounds(1) = 0.05;
    upperBounds(1) = 10;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = -30;
    upperBounds(numel(upperBounds)+1) = 30;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % timeConstantSeconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 2.1;
    lowerBounds(numel(lowerBounds)+1) = 0.0;
    upperBounds(numel(upperBounds)+1) = 5;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    
    % deltaTau (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 30;
    lowerBounds(numel(lowerBounds)+1) = 10;
    upperBounds(numel(upperBounds)+1) = 90;
    paramNames{numel(paramNames)+1} = 'HP delta time constant (msec)';
    
    % leadLagFilterOrder
    initialValues(numel(initialValues)+1) = 4;
    lowerBounds(numel(lowerBounds)+1) = 4;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order';
    

    % lowpass timeConstantSeconds
    initialValues(numel(initialValues)+1) = 2.0;
    lowerBounds(numel(lowerBounds)+1) = 0.2;
    upperBounds(numel(upperBounds)+1) = 5;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';

    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 30;
    paramNames{numel(paramNames)+1} = 'LP filter order';

    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 0.2;
    upperBounds(numel(upperBounds)+1) = 5;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 4;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 30;
    paramNames{numel(paramNames)+1} = 'LP2 filter order';

    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;
    highPassTimeConstantSeconds = theCurrentParams(3)*1e-3;
    deltaTau = theCurrentParams(4)*1e-3;
    lowPassTimeConstantSeconds = theCurrentParams(6)*1e-3;
    lowPass2TimeConstantSeconds = theCurrentParams(8)*1e-3;

    % filter orders must be integer
    theCurrentParams(5) = round(theCurrentParams(5));
    leadLagFilterOrder = theCurrentParams(5);

    theCurrentParams(7) = round(theCurrentParams(7));
    lowPassFilterOrder = theCurrentParams(7);
   
    theCurrentParams(9) = round(theCurrentParams(9));
    lowPass2FilterOrder = theCurrentParams(9);
    
    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    % HighPass (i.e., gain and phase increase with TF) when deltaTau > 0
    theLeadLagFilterTTF = ((1 + 1i * omega * (highPassTimeConstantSeconds+deltaTau)) ./ (1 + 1i * omega * highPassTimeConstantSeconds) ) .^ round(leadLagFilterOrder);

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-lowPassFilterOrder);

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPass2FilterTTF = (1 + 1i * omega * lowPass2TimeConstantSeconds) .^ (-lowPass2FilterOrder);


    theFilterTTF = gain * theDelayFilterTTF .* theLeadLagFilterTTF .* theLowPassFilterTTF .* theLowPass2FilterTTF;
end





function theInnerRetinaTTF = delayHighPassFilterToTransformPhotocurrentTTFtoTargetTTFOLD(...
                temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF, theInnerRetinaTTFfromDirectDeconvolution)


    theTargetTTF = theTargetTTF/max(abs(theTargetTTF));

    useParallel = true;
    multiStartsNum = 32;

    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1100 1000]);
    ax = subplot(1,1,1);

    [~, modelParams.initialValues, ...
        modelParams.lowerBounds, ...
        modelParams.upperBounds, ...
        modelParams.names] = delayHighPassFilter([],temporalFrequencySupportHz);

    modelParams.finalValues = modelParams.initialValues;
    modelParams.scaling = 'linear';

    solverType = 'multi-start';
    solverType = 'global-search';
    %solverType = 'fmincon'

    for iRepeat = 1:3
        figNo = 1000+iRepeat*10;
        objectiveFunctionToMinimize = @(x)theObjectiveFunctionToMinimize(x, temporalFrequencySupportHz, ...
            theTargetTTF, thePhotocurrentsBasedTTF, ax, figNo, modelParams);
       
    
        problem = createOptimProblem('fmincon',...
              'objective', objectiveFunctionToMinimize , ...
              'x0', modelParams.initialValues, ...
              'lb', modelParams.lowerBounds, ...
              'ub', modelParams.upperBounds, ...
              'options', optimoptions(...
                'fmincon',...
                'Display', 'none', ...
                'Algorithm', 'sqp',... % 'sqp', ... % 'interior-point',...
                'GradObj', 'off', ...
                'DerivativeCheck', 'off', ...
                'MaxFunEvals', 10^5, ...
                'MaxIter', 10^4) ...
              );

   

        switch (solverType)
            case 'multi-start'
                % Setup the multi-start solver
                ms = MultiStart(...
                'Display', 'iter', ...
                'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
                'UseParallel', useParallel);
    
                % Run the multi-start
                modelParams.finalValues = run(ms, problem, multiStartsNum);
    
            case 'global-search'
    
                gs = GlobalSearch;
                gs = GlobalSearch(gs,'XTolerance',1e-3,'StartPointsToRun','bounds');
                modelParams.finalValues = run(gs,problem);
    
    
            case 'fmincon'
                modelParams.finalValues = fmincon(problem);
    
            otherwise
                error('Uknown solver type: ''%s''.', solverType)
                
        end
    
        fprintf('Params at iteration: %d\n', iRepeat)
        modelParams.finalValues
        modelParams.initialValues = modelParams.finalValues;

    end % iRepeat



    % Now do it with the photocurrents
    modelParams.initialValues = modelParams.finalValues;
    modelParams.finalValues

    theInnerRetinaTTF = delayHighPassFilter(modelParams.finalValues,temporalFrequencySupportHz);

end



function theResidual = theObjectiveFunctionToMinimizeOLD(theCurrentParams, ...
    temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF,  ax, figNo, modelVariables)

    % 
    [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = delayHighPassFilter(...
        theCurrentParams, temporalFrequencySupportHz);


    modelVariables.finalValues = theCurrentParams;
    RGCMosaicConstructor.visualize.fittedModelParams(ax, modelVariables, 'TTF fit');

    
    performFFTshift = true;
    zeroPaddingLength = 512;
    theCurrentInnerRetinaFilterResponseData = RGCMosaicConstructor.temporalFilterEngine.impulseResponseFunctionFromTTF(...
        theCurrentInnerRetinaFilterTTF, temporalFrequencySupportHz, performFFTshift, zeroPaddingLength);

    hFig = figure(figNo); clf;
    set(hFig, 'Position', [10 10 1000 1000]);
    ax = subplot(2,1,1);
    plot(theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds*1e3, theCurrentInnerRetinaFilterResponseData.amplitude, 'k-');
    set(gca, 'XLim', [500 1500])
    ax = subplot(2,1,2);
    plot(temporalFrequencySupportHz, abs(theTargetTTF), 'r-', 'LineWidth', 1.5);
    hold on
    plot(temporalFrequencySupportHz, abs(thePhotocurrentsBasedTTF), 'k-', 'LineWidth', 1.5);
    plot(temporalFrequencySupportHz, abs(theCurrentInnerRetinaFilterTTF), 'b-', 'LineWidth', 1.5);
    plot(temporalFrequencySupportHz, abs(theCurrentInnerRetinaFilterTTF .* thePhotocurrentsBasedTTF), 'mo', 'LineWidth', 1.5);
    legend({'targetTTF', 'photocurrent TTF', 'inner retina TTF', 'achieved TTF'})
    set(gca, 'XLim', [0.5 200], 'XScale', 'log', 'XTick', [1 3 10 30 100])

    theResidual = norm(theTargetTTF - theCurrentInnerRetinaFilterTTF.*thePhotocurrentsBasedTTF);
end


function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayHighPassFilter(theCurrentParams, temporalFrequencySupportHz)
    
    % gain
    initialValues(1) = 5;
    lowerBounds(1) = 1;
    upperBounds(1) = 5;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = -40;
    lowerBounds(numel(lowerBounds)+1) = -60;
    upperBounds(numel(upperBounds)+1) = 60;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % highpass gain
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 5.0;
    paramNames{numel(paramNames)+1} = 'highpaass gain';
  
    % highpass timeConstant msec
    initialValues(numel(initialValues)+1) = 0.2;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 0.3;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    

    % lowpass time cosntat msec
    initialValues(numel(initialValues)+1) = 0.1;
    lowerBounds(numel(lowerBounds)+1) = 0.0;
    upperBounds(numel(upperBounds)+1) = 2;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';
    
    % nLtL
    initialValues(numel(initialValues)+1) = 50;
    lowerBounds(numel(lowerBounds)+1) = 45;
    upperBounds(numel(upperBounds)+1) = 55;
    paramNames{numel(paramNames)+1} = 'nLtL';
    

    % high pass FilterOrder
    initialValues(numel(initialValues)+1) = 50;
    lowerBounds(numel(lowerBounds)+1) = 30;
    upperBounds(numel(upperBounds)+1) = 200;
    paramNames{numel(paramNames)+1} = 'HP filter order';



    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    % Get params
    gain = theCurrentParams(1);                              % Responsitivity at 0 Hz
    conductionDelaySeconds = theCurrentParams(2)*1e-3;       % D
    highPassGain = theCurrentParams(3);                      % Hs
    highPassTimeConstantSeconds = theCurrentParams(4)*1e-3;  % Tau_s (Benardete and Kaplan (1992a) varied this for different contrast levels)
    lowPassTimeConstantSeconds = theCurrentParams(5)*1e-3;   % Tau_l
    
    nLowPassStagesNum = round(theCurrentParams(6)/theCurrentParams(5)) ;               % NlTl
    nHighPassStagesNum = theCurrentParams(7);                % always 1 in Benardete & Kaplan (1992a)

    % filter orders must be integer
    theCurrentParams(7) = max([1 round(theCurrentParams(7))]);


    % Circular frequency in radians
    omega = 2 * pi * temporalFrequencySupportHz;

    % Delay filter
    theDelayFilterTTF = exp(-1i * omega * conductionDelaySeconds);

    % 1-stage high-pass filter
    theHighPassFilterTTF = 1 - highPassGain * (1 + 1i * omega * highPassTimeConstantSeconds) .^ (-nHighPassStagesNum);

    % N-stage low-pass filter
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-nLowPassStagesNum);

    theFilterTTF = gain * theDelayFilterTTF .* theHighPassFilterTTF .* theLowPassFilterTTF;


end


