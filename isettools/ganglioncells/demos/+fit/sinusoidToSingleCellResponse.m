function [responseTimeAxisHR, fittedSinusoid, fittedParams] = sinusoidToSingleCellResponse(...
            responseTimeAxis, cellResponse, temporalFrequencyHz)
        
        sinFunction = @(params,time)(params(1) * sin(2.0*pi*temporalFrequencyHz*time - params(2)));
        opts.RobustWgtFun = []; %'talwar';
        opts.MaxIter = 1000;
        
        % Subtract mean response
        meanResponse = mean(cellResponse);
        cellResponse = cellResponse - meanResponse;

        % Fit sinusoid
        initialParams(1) = prctile(cellResponse, 85); % gain
        initialParams(2) = 0.0;                       % phase
        fittedParams = nlinfit(responseTimeAxis, cellResponse, sinFunction, initialParams, opts);
        if (fittedParams(1) < 0)
            fittedParams(1) = -fittedParams(1);
            fittedParams(2) = fittedParams(2) + pi;
        end
            
        % Compute high resolution fitted response
        responseTimeAxisHR = linspace(responseTimeAxis(1), responseTimeAxis(end), 100);
        fittedSinusoid = sinFunction(fittedParams, responseTimeAxisHR) + meanResponse;
        
        % Return gain and phase of fitted sinusoid
        fittedParams(1) = fittedParams(1);
        fittedParams(2) = fittedParams(2)/pi*180;
end

