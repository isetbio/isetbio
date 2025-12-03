function t_temporalImpulseResponse

	% Compute TTF
	temporalFrequencySupportHz = linspace(0.1, 64, 16);

    nL_tL_product = 48;

    % Center TTF params (ON-center cell data from Figure 6)
    centerParams(1) = 184.2;        % gain (A)

    centerParams(2) = 0.69;         % high-pass gain (Hs)
    centerParams(3) = 18.61;        % high-pass time constant (msec) (Ts)

    centerParams(4) = 38;           % low-pass stages num (Nl)
    centerParams(5) = nL_tL_product/centerParams(4);         % low-pass time constant (msec) (Tl)

    centerParams(6) = 4.0;          % delay (msec) (D)


    % Surround TTF params (ON-center cell data from Figure 6)
    surroundParams(1) = 125.33;       % gain (A)
    surroundParams(2) = 0.56;         % high-pass gain (Hs)
    surroundParams(3) = 33.28;        % high-pass time constant (msec) (Ts)
    surroundParams(4) = 124;          % low-pass stages num (Nl)
    surroundParams(5) = nL_tL_product/surroundParams(4);         % low-pass time constant (msec) (Tl)
    surroundParams(6) = 4.03;         % delay (msec) (D)


    if (1==2)
        % Values from Figure 7
        centerParams(1) = 184.2;        % gain (A)
        centerParams(2) = 0.70;         % high-pass gain (Hs)
        centerParams(3) = 37.30;        % high-pass time constant (msec) (Ts)
        centerParams(4) = 32;           % low-pass stages num (Nl)
        centerParams(5) = nL_tL_product/centerParams(4);         % low-pass time constant (msec) (Tl)
        centerParams(6) = 4.0;          % delay (msec) (D)

        % Values from Figure 7
        surroundParams(1) = centerParams(1)*0.64;       % gain (A)
        surroundParams(2) = 0.41;         % high-pass gain (Hs)
        surroundParams(3) = 42.49;        % high-pass time constant (msec) (Ts)
        surroundParams(4) = 93;          % low-pass stages num (Nl)
        surroundParams(5) = nL_tL_product/surroundParams(4);         % low-pass time constant (msec) (Tl)
        surroundParams(6) = 4.0;          % delay (msec) (D)
    end


    % Compute temporal transfer functions
	theCenterTTF = highPassNstageLowPassTTF(centerParams, temporalFrequencySupportHz);
    theSurroundTTF = highPassNstageLowPassTTF(surroundParams, temporalFrequencySupportHz);


    % Compute the product TTF (just for testing) cascade-filter
    a = real(theCenterTTF);
    b = imag(theCenterTTF);
    c = real(theSurroundTTF);
    d = imag(theSurroundTTF);
    theProductTTF = (a.*c - b.*d) + 1i * (a.*d + b.*c);
    theProductTTF = theProductTTF ./ norm(theCenterTTF);

    gain = struct(...
        'initial', centerParams(1), ...
        'low', centerParams(1)/10, ...
        'high', centerParams(1)*10 ...
        );

    highPassGain = struct(...
        'initial',centerParams(2), ...
        'low', 0, ...
        'high', 1 ...
        );

    highPassTimeConstant = struct(...
        'initial',centerParams(3), ...
        'low', 0, ...
        'high', 100 ...
        );

    lowPassStagesNum = struct(...
        'initial',centerParams(4), ...
        'low', 1, ...
        'high', 300 ...
        );

    lowPassTimeConstant = struct(...
        'initial',centerParams(5), ...
        'low', 0.1, ...
        'high', 10 ...
        );

    delay = struct(...
        'initial', centerParams(6), ...
        'low', 0, ...
        'high', 400 ...
        );

    TTFparams.initialValues = [gain.initial   highPassGain.initial    highPassTimeConstant.initial    lowPassStagesNum.initial  lowPassTimeConstant.initial      delay.initial];
    TTFparams.lowerBounds   = [gain.low       highPassGain.low        highPassTimeConstant.low        lowPassStagesNum.low      lowPassTimeConstant.low          delay.low];
    TTFparams.upperBounds   = [gain.high      highPassGain.high       highPassTimeConstant.high       lowPassStagesNum.high     lowPassTimeConstant.high         delay.high];
    TTFparams.names         = {'gain',        'highPassGain',         'high pass tau (msec)',  'low pass stages',       'low pass tau (msec)',    'delay (msec)'};
    TTFparams.scaling       = {'log',         'linear',               'linear',                       'linear',                 'linear',                        'linear'};

    hFig = figure(100);
    axModelParams = subplot(1,1,1);

    % Fit theProductTTF
    [TTFparams, theFitttedToTheProductTTF] = highPassNstageLowPassTTFtoArbitraryTTF(...
        theProductTTF, temporalFrequencySupportHz, TTFparams, axModelParams);


    % Convert to temporal impulse response functions
    theCenterImpulseResponseData = ...
        temporalTransferFunctionToImpulseResponseFunction(theCenterTTF, temporalFrequencySupportHz);
    theSurroundImpulseResponseData  = ...
        temporalTransferFunctionToImpulseResponseFunction(theSurroundTTF, temporalFrequencySupportHz);

    % Compute the IR of the cascade filter
    theProductImpulseResponseData = ...
        temporalTransferFunctionToImpulseResponseFunction(theProductTTF, temporalFrequencySupportHz);

    theFittedToTheProductImpulseResponseData = ...
        temporalTransferFunctionToImpulseResponseFunction(theFitttedToTheProductTTF, temporalFrequencySupportHz);


    % Plot temporal transfer functions
    figure(1); clf;
    subplot(1,2,1)
    p1 = plot(temporalFrequencySupportHz, abs(theCenterTTF), 'ro-');
    hold on;
    p2 = plot(temporalFrequencySupportHz, abs(theSurroundTTF), 'bo-');
    p3 = plot(temporalFrequencySupportHz, abs(theProductTTF), 'mo-');
    p4 = plot(temporalFrequencySupportHz, abs(theFitttedToTheProductTTF), 'k.');
    set(gca, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.25 0.5 1 2 4 8 16 32 64 128]);
    legend([p1 p2 p3 p4], {'center', 'surround', 'center x surround', 'fitted to center x surround'});
    grid on

    subplot(1,2,2)
    plot(temporalFrequencySupportHz, unwrap(angle(theCenterTTF))/pi*180, 'ro-');
    hold on
    plot(temporalFrequencySupportHz, unwrap(angle(theCenterTTF)-angle(theSurroundTTF))/pi*180, 'ko--');
    set(gca, 'XScale', 'log', 'XLim', [0.25 32], 'XTick', [0.25 0.5 1 2 4 8 16 32], 'YLim', [0 180], 'YTick', -360:30:360);
    ylabel('center-surround phase difference (degs)');
    grid on

    % Plot impulse response functions
    figure(2); clf;
    p1 = plot(theCenterImpulseResponseData.temporalSupportSeconds*1e3, theCenterImpulseResponseData.weights, 'ro-');
    hold on;
    p2 = plot(theSurroundImpulseResponseData.temporalSupportSeconds*1e3, theSurroundImpulseResponseData.weights, 'bo-');
    p3 = plot(theProductImpulseResponseData.temporalSupportSeconds*1e3, theProductImpulseResponseData.weights, 'mo-');
    p4 = plot(theFittedToTheProductImpulseResponseData.temporalSupportSeconds*1e3, theFittedToTheProductImpulseResponseData.weights, 'k.-');
    legend([p1 p2 p3 p4], {'center', 'surround', 'center+surround', 'fitted to center x surround'});
    xlabel('time (msec)')
    set(gca, 'XLim', [0 200], 'XTick', 0:10:1000);
    grid on

end


% Fit the highPassNstageLowPassTTF to an arbitrary TTF
function  [TTFparams, theFittedTTF] = highPassNstageLowPassTTFtoArbitraryTTF(theArbitraryComplexTTF, temporalFrequencySupportHz, TTFparams, ax)

    objective = @(x)highPassNstageLowPassTTFresidual(x, theArbitraryComplexTTF, temporalFrequencySupportHz, ax, TTFparams);

    % Multi-start
    problem = createOptimProblem('fmincon',...
          'objective', objective, ...
          'x0', TTFparams.initialValues, ...
          'lb', TTFparams.lowerBounds, ...
          'ub', TTFparams.upperBounds, ...
          'options', optimoptions(...
            'fmincon',...
            'Display', 'none', ...
            'Algorithm', 'interior-point',... % 'sqp', ... % 'interior-point',...
            'GradObj', 'off', ...
            'DerivativeCheck', 'off', ...
            'MaxFunEvals', 10^5, ...
            'MaxIter', 10^4) ...
          );

    ms = MultiStart(...
          'Display', 'iter', ...
          'StartPointsToRun','bounds-ineqs', ...  % run only initial points that are feasible with respect to bounds and inequality constraints.
          'UseParallel', true);

    multiStartsNum = 32;

    % Run the multi-start
    TTFparams.finalValues = run(ms, problem, multiStartsNum)


    theFittedTTF = highPassNstageLowPassTTF(TTFparams.finalValues, temporalFrequencySupportHz);

end




function theImpulseResponseFunctionData = ...
        temporalTransferFunctionToImpulseResponseFunction(theTTF, temporalFrequencySupportHz)
    % Convert single-sided spectrum to double sided
    theDoubleSidedTTF = [theTTF(1) theTTF(2:end)/2 fliplr(conj(theTTF(2:end)))/2];

    theImpulseResponse = ifft(theDoubleSidedTTF, 'symmetric');

    fMax = max(temporalFrequencySupportHz);
    dtSeconds = 1/(2*fMax)
    theTemporalSupportSeconds = (1:numel(theImpulseResponse)) * dtSeconds;

    theImpulseResponseFunctionData.weights = theImpulseResponse;
    theImpulseResponseFunctionData.temporalSupportSeconds = theTemporalSupportSeconds;
end


function theTTF = highPassNstageLowPassTTF(params, temporalFrequencySupportHz)

	% Get params
    gain = params(1);
    highPassGain = params(2);
    highPassTimeConstantSeconds = params(3)*1e-3;
    lowPassStagesNum = round(params(4));
    lowPassTimeConstantSeconds = params(5)*1e-3;
    delaySeconds = params(6)*1e-3;
    omega = temporalFrequencySupportHz * (2 * pi);

    % The TwoStageTTF model in the frequency domain
    theTTF = ...
    	gain * exp(-1i * omega * delaySeconds) .* ...
    	(1 - highPassGain * 1./(1 + 1i * omega * highPassTimeConstantSeconds)) .* ...
    	(1./(1 + 1i * omega * lowPassTimeConstantSeconds)) .^ (lowPassStagesNum);

end


function theResidual = highPassNstageLowPassTTFresidual(theCurrentParams, theTTFtoFit, temporalFrequencySupportHz, ax, modelVariables)

    theResidual = norm(highPassNstageLowPassTTF(theCurrentParams, temporalFrequencySupportHz) - theTTFtoFit);

    modelVariables.finalValues = theCurrentParams;
    RGCMosaicConstructor.visualize.fittedModelParams(ax, modelVariables, 'TTF fit');

end
