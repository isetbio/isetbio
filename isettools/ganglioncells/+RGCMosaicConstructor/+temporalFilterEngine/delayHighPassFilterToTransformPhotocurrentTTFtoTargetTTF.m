%
% RGCMosaicConstructor.temporalFilterEngine.delayHighPassFilterToTransformPhotocurrentTTFtoTargetTTF
%

function theInnerRetinaTTF = delayHighPassFilterToTransformPhotocurrentTTFtoTargetTTF(...
                temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF, ...
                filterType, solverType, multiStartsNum, useParallel)


    theTargetTTF = theTargetTTF / max(abs(theTargetTTF));

    hFig = figure(1000); clf;
    set(hFig, 'Position', [10 10 1100 1000]);
    ax = subplot(1,1,1);

    switch (filterType)
        case 'delayLeadLagFilter'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = delayLeadLagFilter([],temporalFrequencySupportHz);

        case 'delayLeadLagFilter2'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = delayLeadLagFilter2([],temporalFrequencySupportHz);

        case 'delayHighPassFilter'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = delayHighPassFilter([],temporalFrequencySupportHz);

        case 'asymmetricBandPassFilter'
            [~, modelParams.initialValues, ...
                modelParams.lowerBounds, ...
                modelParams.upperBounds, ...
                modelParams.names] = asymmetricBandPassFilter([],temporalFrequencySupportHz);
        otherwise
            error('Unknown filterTye: ''%s''.', filterType);
    end

    modelParams.finalValues = modelParams.initialValues;
    modelParams.scaling = 'linear';


    showCurrentParamValuesPlot = true;
    showFitResult = true;

    hFig = figure(2000); clf;
    set(hFig, 'Position', [500 300 1100 1000]);
    ax = subplot(1,1,1);
    figNo = 2001;

    

    objectiveFunctionToMinimize = @(x)theObjectiveFunctionToMinimize(x, temporalFrequencySupportHz, ...
        theTargetTTF, thePhotocurrentsBasedTTF, filterType, ax, figNo, modelParams, ...
        showCurrentParamValuesPlot, showFitResult);

    
    residualsSequence = [];
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
                modelParams.finalValues = run(ms, problem, multiStartsNum)
    
            case 'global-search'
    
                gs = GlobalSearch;
                gs = GlobalSearch(gs,'XTolerance',1e-3,'StartPointsToRun','bounds');
                modelParams.finalValues = run(gs,problem)
    
    
            case 'fmincon'
                modelParams.finalValues = fmincon(problem)
    
            otherwise
                error('Uknown solver type: ''%s''.', solverType)
                
     end % switch (solver)

  

    % For lead-lag filter order fixed to 1.0, modelParams.finalValues were found as follows:
    % 0.1306 -7.6329 1.1930 150 (MAX) 1 0.9796 20.9183 0.9796 13.9997

    % For lead-lag filter order fixed to 2.0, modelParams.finalValues were found as follows:
    % 0.0627 -0.8099 1.2076 71.8375   2 1.2076 14.4942 1.2076 19.8082

    % For lead-lag filter order fixed to 3.0, modelParams.finalValues were found as follows:
    %0.1335  -9.1854 1.6700 29.5640   3  1.6700 4.1152 1.6700 28.7999

    % For lead-lag filter order fixed to 4.0, modelParams.finalValues were found as follows:
    % 0.1853 -29.4469 1.6509 20.0019  4 1 .6509 18.0088 1.6508 31.3295


    % For lead-lag filter order fixed to 5.0, modelParams.finalValues were found as follows:
    % 0.1942 -5.9163 2.8167 15.7018 5 2.8167 8.9332 2.8166 11.6619

    % For lead-lag filter order fixed to 6.0, modelParams.finalValues were found as follows:
    % 0.2153 -13.8488 2.6219 13.5265 6 2.6218 24.4346 5 1.6532

    % For lead-lag filter order fixed to 4.0 (256 multi-starts), modelParams.finalValues were found as follows:
    % 0.1738 -9.5041 2.2091 20.0430 4 2.2091 23.6671 2.2091 3.2103

    switch (filterType)
        case 'delayLeadLagFilter'
            theInnerRetinaTTF = delayLeadLagFilter(modelParams.finalValues,temporalFrequencySupportHz);

        case 'delayLeadLagFilter2'
            theInnerRetinaTTF = delayLeadLagFilter2(modelParams.finalValues,temporalFrequencySupportHz);

        case 'delayHighPassFilter'
            theInnerRetinaTTF = delayHighPassFilter(modelParams.finalValues,temporalFrequencySupportHz);

        case 'asymmetricBandPassFilter'
            theInnerRetinaTTF = asymmetricBandPassFilter(modelParams.finalValues,temporalFrequencySupportHz);
    end


    % Nested function
    function theResidual = theObjectiveFunctionToMinimize(theCurrentParams, ...
        temporalFrequencySupportHz, theTargetTTF, thePhotocurrentsBasedTTF,  ...
        filterType, ax, figNo, modelParams, showCurrentParamValuesPlot, showFitResults)

        switch (filterType)
            case 'delayLeadLagFilter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    delayLeadLagFilter(theCurrentParams, temporalFrequencySupportHz);
    
            case 'delayLeadLagFilter2'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    delayLeadLagFilter2(theCurrentParams, temporalFrequencySupportHz);
    
            case 'delayHighPassFilter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    delayHighPassFilter(theCurrentParams, temporalFrequencySupportHz);
    
            case 'asymmetricBandPassFilter'
                [theCurrentInnerRetinaFilterTTF, ~, ~, ~, ~, theCurrentParams] = ...
                    asymmetricBandPassFilter(theCurrentParams, temporalFrequencySupportHz);
        end

        theResidual = norm(theTargetTTF - theCurrentInnerRetinaFilterTTF.*thePhotocurrentsBasedTTF);
    
        residualsSequence(numel(residualsSequence)+1) = theResidual;
        if (theResidual == min(residualsSequence(:)))
            figure(55)
            plot(1:numel(residualsSequence),residualsSequence, 'ks-' )
            drawnow;
            save(sprintf('OptimizationAtIteration_%d.mat', numel(residualsSequence)), ...
                'theCurrentParams', ...
                'temporalFrequencySupportHz', ...
                'theCurrentInnerRetinaFilterTTF', ...
                'thePhotocurrentsBasedTTF', ...
                'theTargetTTF')
        end

    
        if (showCurrentParamValuesPlot)
            % Show the current param values plot
            modelParams.finalValues = theCurrentParams;
            RGCMosaicConstructor.visualize.fittedModelParams(ax, modelParams, 'TTF fit');
        end
    
    
        if (showFitResults)
            theCurrentInnerRetinaFilterResponseData = RGCMosaicConstructor.temporalFilterEngine.sampledTTFtoTemporalImpulseFunction(...
                theCurrentInnerRetinaFilterTTF, temporalFrequencySupportHz);
    
            hFig = figure(figNo); clf;
            set(hFig, 'Position', [10 10 1000 1250]);
            ax = subplot(1,2,1);
            plot(theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds*1e3, theCurrentInnerRetinaFilterResponseData.amplitude, 'k-');
            m = numel(theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds);
            set(gca, 'XLim', theCurrentInnerRetinaFilterResponseData.temporalSupportSeconds(m)+[-100 100]);
            ax = subplot(1,2,2);
            plot(temporalFrequencySupportHz, abs(theTargetTTF), 'r-');
            hold on
            plot(temporalFrequencySupportHz, abs(thePhotocurrentsBasedTTF), 'k-');
            plot(temporalFrequencySupportHz, abs(theCurrentInnerRetinaFilterTTF), 'b-');
            plot(temporalFrequencySupportHz, abs(theCurrentInnerRetinaFilterTTF .* thePhotocurrentsBasedTTF), 'm--');
            set(gca, 'XLim', [0.5 200], 'YLim', [0 2*max(abs(theTargetTTF))], 'XScale', 'log', 'XTick', [1 3 10 30 100])
        end
    end % Nested function

end


function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    asymmetricBandPassFilter(theCurrentParams, temporalFrequencySupportHz)

    % gain
    initialValues(1) = 0.13;
    lowerBounds(1) = 0.5;
    upperBounds(1) = 30;
    paramNames{1} = 'gain';

    % delaySeconds
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = -60;
    upperBounds(numel(upperBounds)+1) = 60;
    paramNames{numel(paramNames)+1} = 'delay (msec)';

    % tau lowpass
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 3;
    paramNames{numel(paramNames)+1} = 'tau lowpass (msec)';


    % tau Low / tau High difference (must be > 1)
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 100;
    paramNames{numel(paramNames)+1} = 'tau lowpass / tau highpass(msec)';


    % order (lowpass) (must be >= 0)
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 2;
    upperBounds(numel(upperBounds)+1) = 2;
    paramNames{numel(paramNames)+1} = 'order lowpass';

    % order (highpass) (must be >= 1)
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'order highpass';

    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end

 
    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;

    tauLowPassSeconds = theCurrentParams(3)*1e-3;
    tauHighPassSeconds = max([1 theCurrentParams(3)*theCurrentParams(4)])*1e-3;
    
    % WE MUST HAVE: tauHighPassSeconds < tauLowPassSeconds,
    % so that highPassFrequency > lowPassFrequency

    orderLowPass = max([0 round(theCurrentParams(5))]);
    orderHighPass = max([1 round(theCurrentParams(6))]);

    
    highPassFrequency = 1/tauHighPassSeconds;
    lowPassFrequency = 1/tauLowPassSeconds;
    centerFrequency = sqrt(highPassFrequency * lowPassFrequency);


    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    theLowPassFilterTTF = (1/tauLowPassSeconds ./ (1i * omega + 1/tauLowPassSeconds)) .^ orderLowPass;
    theHighPassFilterTTF = ((1i * omega) ./ (1i * omega + 1/tauHighPassSeconds)) .^ orderHighPass;
   

    theFilterTTF = gain * theDelayFilterTTF .* theHighPassFilterTTF .* theLowPassFilterTTF;

end

function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayLeadLagFilter(theCurrentParams, temporalFrequencySupportHz)
    

    % gain
    initialValues(1) = 0.13;
    lowerBounds(1) = 0.05;
    upperBounds(1) = 10;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = -30;
    lowerBounds(numel(lowerBounds)+1) = -60;
    upperBounds(numel(upperBounds)+1) = 60;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % timeConstantSeconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 0.2;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    
    % deltaTau (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 20;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 200;
    paramNames{numel(paramNames)+1} = 'HP delta time constant (msec)';
    
    % leadLagFilterOrder
    initialValues(numel(initialValues)+1) = 3;
    lowerBounds(numel(lowerBounds)+1) = 2;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order';
    

    % lowpass timeConstantSeconds
    initialValues(numel(initialValues)+1) = 0.65;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';

    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 5;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 50;
    paramNames{numel(paramNames)+1} = 'LP filter order';

    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 100;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';

    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 10;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 100;
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
    %theCurrentParams(5) = round(theCurrentParams(5));
    leadLagFilterOrder = theCurrentParams(5);

    %theCurrentParams(7) = round(theCurrentParams(7));
    lowPassFilterOrder = theCurrentParams(7);
   
    %theCurrentParams(9) = round(theCurrentParams(9));
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


function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayLeadLagFilter2(theCurrentParams, temporalFrequencySupportHz)
    
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
    
    % timeConstant1Seconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 2.1;
    lowerBounds(numel(lowerBounds)+1) = 0.0;
    upperBounds(numel(upperBounds)+1) = 25;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    
    % timeConstant2Seconds (tau2 in Purpura, Tranchina, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 2.1;
    lowerBounds(numel(lowerBounds)+1) = 0.0;
    upperBounds(numel(upperBounds)+1) = 25;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';

    % deltaTau1 (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 60;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 300;
    paramNames{numel(paramNames)+1} = 'HP delta time constant 1 (msec)';
    
    % deltaTau2 (timeConstantSeconds + deltaTau = tau1 in Purpura, Tranchine, Kaplan Shapley 1990)
    initialValues(numel(initialValues)+1) = 60;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 300;
    paramNames{numel(paramNames)+1} = 'HP delta time constant 2 (msec)';

    % leadLagFilterOrder1
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order 1';
    
    % leadLagFilterOrder2
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'lead-lag filter order 2';



    % lowpass timeConstantSeconds
    initialValues(numel(initialValues)+1) = 2.0;
    lowerBounds(numel(lowerBounds)+1) = 0.2;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';

    % lowpass2 timeConstantSeconds
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 20;
    paramNames{numel(paramNames)+1} = 'LP2 time constant (msec)';

    % lowpass FilterOrder
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 2;
    paramNames{numel(paramNames)+1} = 'LP filter order';

    
    % lowpass2 FilterOrder
    initialValues(numel(initialValues)+1) = 2;
    lowerBounds(numel(lowerBounds)+1) = 2;
    upperBounds(numel(upperBounds)+1) = 4;
    paramNames{numel(paramNames)+1} = 'LP2 filter order';

    if (isempty(theCurrentParams))
        theFilterTTF = [];
        return;
    end


    gain = theCurrentParams(1);
    delaySeconds = theCurrentParams(2)*1e-3;
    highPassTimeConstant1Seconds = theCurrentParams(3)*1e-3;
    highPassTimeConstant2Seconds = theCurrentParams(4)*1e-3;
    deltaTau1 = theCurrentParams(5)*1e-3;
    deltaTau2 = theCurrentParams(6)*1e-3;

    % filter orders must be integer
    theCurrentParams(7) = round(theCurrentParams(7));
    leadLagFilter1Order = theCurrentParams(7);

    % filter orders must be integer
    theCurrentParams(8) = round(theCurrentParams(8));
    leadLagFilter2Order = theCurrentParams(8);


    lowPassTimeConstantSeconds = theCurrentParams(9)*1e-3;
    lowPass2TimeConstantSeconds = theCurrentParams(10)*1e-3;

    
    % filter orders must be integer
    theCurrentParams(11) = round(theCurrentParams(11));
    lowPassFilterOrder = theCurrentParams(11);
   
    % filter orders must be integer
    theCurrentParams(12) = round(theCurrentParams(12));
    lowPass2FilterOrder = theCurrentParams(12);
    
    omega = 2 * pi * temporalFrequencySupportHz;
    theDelayFilterTTF = exp(-1i * omega * delaySeconds);

    % HighPass (i.e., gain and phase increase with TF) when deltaTau > 0
    theLeadLagFilter1TTF = ((1 + 1i * omega * (highPassTimeConstant1Seconds+deltaTau1)) ./ (1 + 1i * omega * highPassTimeConstant1Seconds) ) .^ round(leadLagFilter1Order);

    % HighPass (i.e., gain and phase increase with TF) when deltaTau > 0
    theLeadLagFilter2TTF = ((1 + 1i * omega * (highPassTimeConstant2Seconds+deltaTau2)) ./ (1 + 1i * omega * highPassTimeConstant2Seconds) ) .^ round(leadLagFilter2Order);


    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-lowPassFilterOrder);

    % N-stage low-pass filter with lowPassTimeConstantSeconds time constant 
    theLowPass2FilterTTF = (1 + 1i * omega * lowPass2TimeConstantSeconds) .^ (-lowPass2FilterOrder);


    theFilterTTF = gain * theDelayFilterTTF .* theLeadLagFilter1TTF .* theLeadLagFilter2TTF .* theLowPassFilterTTF .* theLowPass2FilterTTF;
end





function [theFilterTTF, initialValues, lowerBounds, upperBounds, paramNames, theCurrentParams] = ...
    delayHighPassFilter(theCurrentParams, temporalFrequencySupportHz)
    
    % gain
    initialValues(1) = 5;
    lowerBounds(1) = 0.1;
    upperBounds(1) = 20;
    paramNames{1} = 'gain';
    
    % delaySeconds
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = -0;
    upperBounds(numel(upperBounds)+1) = 0;
    paramNames{numel(paramNames)+1} = 'delay (msec)';
    
    % highpass gain
    initialValues(numel(initialValues)+1) = 0.0;
    lowerBounds(numel(lowerBounds)+1) = 0.0;
    upperBounds(numel(upperBounds)+1) = 0.0;
    paramNames{numel(paramNames)+1} = 'highpaass gain';
  
    % highpass timeConstant msec
    initialValues(numel(initialValues)+1) = 0.1;
    lowerBounds(numel(lowerBounds)+1) = 0.1;
    upperBounds(numel(upperBounds)+1) = 0.1;
    paramNames{numel(paramNames)+1} = 'HP time constant (msec)';
    

    % lowpass time cosntat msec
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 0;
    paramNames{numel(paramNames)+1} = 'LP time constant (msec)';
    
    % lowpass filter order
    initialValues(numel(initialValues)+1) = 1;
    lowerBounds(numel(lowerBounds)+1) = 1;
    upperBounds(numel(upperBounds)+1) = 1;
    paramNames{numel(paramNames)+1} = 'LP filter order';
    

    % high pass FilterOrder
    initialValues(numel(initialValues)+1) = 0;
    lowerBounds(numel(lowerBounds)+1) = 0;
    upperBounds(numel(upperBounds)+1) = 0;
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
    
    
    % filter orders must be integer
    theCurrentParams(6) = max([1 round(theCurrentParams(6))]);
    theCurrentParams(7) = max([1 round(theCurrentParams(7))]);

    nLowPassStagesNum = theCurrentParams(6) ;               % NlTl
    nHighPassStagesNum = theCurrentParams(7);                % always 1 in Benardete & Kaplan (1992a)


    % Circular frequency in radians
    omega = 2 * pi * temporalFrequencySupportHz;

    % Delay filter
    theDelayFilterTTF = 1+0*exp(-1i * omega * conductionDelaySeconds);

    % 1-stage high-pass filter
    theHighPassFilterTTF = 1 - highPassGain * (1 + 1i * omega * highPassTimeConstantSeconds) .^ (-nHighPassStagesNum);

    % N-stage low-pass filter
    theLowPassFilterTTF = (1 + 1i * omega * lowPassTimeConstantSeconds) .^ (-nLowPassStagesNum);

    theFilterTTF = 1+0*gain * theDelayFilterTTF .* theHighPassFilterTTF .* theLowPassFilterTTF;


end


