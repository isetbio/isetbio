
function [responseAmplitude, responsePhase, responseTimeAxisHR, fittedResponses] = ...
    fitSinusoidalModulationToResponseTimeCourse(responses, responsesStDev, responseTimeAxis, ...
        stimulusTemporalFrequencyHz, visualizeFits, spatialFrequency, maxSpikeRate)
        % Compute response modulations
        rgcsNum = size(responses,1);
        responseAmplitude = zeros(rgcsNum,1);
        responsePhase = zeros(rgcsNum,1);
        responseTimeAxisHR = responseTimeAxis(1):0.5/1000.0:responseTimeAxis(end);
        fittedResponses = zeros(rgcsNum, numel(responseTimeAxisHR));
        
        
        if (visualizeFits)
            hFig = figure(333);
            set(hFig, 'Position', [10 10 1400 800]);
        end
        
        sinFunction = @(params,time)(params(1) * sin(2.0*pi*stimulusTemporalFrequencyHz*time - params(2)));
        opts.RobustWgtFun = []; %'talwar';
        opts.MaxIter = 1000;
            
        for iRGC = 1:rgcsNum
            response = squeeze(responses(iRGC,:));
            weights = 1.0 ./ squeeze(responsesStDev(iRGC,:));
            % Fit sinusoid with phase, amplitude free params
            prctiles = prctile(response,[15 85]);
            initialParams(1) = prctiles(2); % amplitude
            initialParams(2) = 0.0; % phase
            fittedParams = nlinfit(responseTimeAxis,response, sinFunction,initialParams, opts, 'Weights', weights);
            if (fittedParams(1)<0)
                fittedParams(1) = -fittedParams(1);
                fittedParams(2) = fittedParams(2) + pi;
            end
            responseAmplitude(iRGC) = fittedParams(1);     % amplitude
            responsePhase(iRGC) = fittedParams(2)/pi*180;   % phase
            responseHR = sinFunction(fittedParams,responseTimeAxisHR);
            fittedResponses(iRGC,:) = responseHR;
            if (visualizeFits)
                errorbar(responseTimeAxis, response, 1./weights, 'ko-', 'MarkerFaceColor', [0.8 0.8 0.8]); hold on;
                plot(responseTimeAxisHR, responseHR, 'r-', 'LineWidth', 1.5); hold off;
                set(gca, 'YLim', maxSpikeRate*[-1 1], 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]);
                title(sprintf('IRGC: %d (sf = %2.2f c/deg)', iRGC, spatialFrequency))
                drawnow;
            end
        end
        pause(1.0)
end
