function [responseAmplitude, responseAmplitudeSE, responsePhase, responseTimeAxisHR, fittedResponses] = ...
    fitSinusoidalModulationToResponseTimeCourse(responses, responsesStDev, responseTimeAxis, ...
        stimulusTemporalFrequencyHz, spatialFrequency, maxSpikeRate, ...
        visualizeIndividualFits, targetRGC, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, exportDir)
    
        % Compute response modulations
        rgcsNum = size(responses,1);
        responseAmplitude = zeros(rgcsNum,1);
        responseAmplitudeSE = zeros(rgcsNum,1);
        responsePhase = zeros(rgcsNum,1);
        responseTimeAxisHR = responseTimeAxis(1):0.5/1000.0:responseTimeAxis(end);
        fittedResponses = zeros(rgcsNum, numel(responseTimeAxisHR));
        
        sinFunction = @(params,time)(params(3) + params(1) * sin(2.0*pi*stimulusTemporalFrequencyHz*time - params(2)));
        opts.RobustWgtFun = []; %'talwar';
        opts.MaxIter = 1000;
            
        if (visualizeIndividualFits)
            plotlabOBJ = setupPlotLab(0);
        end
    
        for iRGC = 1:rgcsNum
            response = squeeze(responses(iRGC,:));
            meanResponse = mean(response);
            weights = 1.0 ./ squeeze(responsesStDev(iRGC,:));
            % Fit sinusoid with phase, amplitude free params
            prctiles = prctile(response-meanResponse,[15 85]);
            initialParams(1) = prctiles(2); % amplitude
            initialParams(2) = 0.0; % phase
            initialParams(3) = mean(response);  % mean
            [fittedParams,~,~,varCovarianceMatrix,~] = nlinfit(responseTimeAxis,response, sinFunction,initialParams, opts, 'Weights', weights);
            if (fittedParams(1)<0)
                fittedParams(1) = -fittedParams(1);
                fittedParams(2) = fittedParams(2) + pi;
            end
            % standard error of the mean
            fittedParamsSE = sqrt(diag(varCovarianceMatrix));
            responseAmplitude(iRGC) = fittedParams(1);     % amplitude
            responsePhase(iRGC) = fittedParams(2)/pi*180;   % phase
            responseAmplitudeSE(iRGC) = fittedParamsSE(1); % SE of the amplitude
            responseHR = sinFunction(fittedParams,responseTimeAxisHR);
            fittedResponses(iRGC,:) = responseHR;
            
            if (visualizeIndividualFits) && ((iRGC == targetRGC) || (isempty(targetRGC)))
                hFig = figure(333); clf;
                errorbar(responseTimeAxis, response, 1./weights, 'ko-', 'LineWidth', 1.0, 'Color', [0.4 0.4 0.4]); hold on;
                line(responseTimeAxisHR, responseHR, 'Color', [0 1 1], 'LineWidth', 3); hold off;
                set(gca, 'YLim', maxSpikeRate*[0 1], 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]);
                title(sprintf('RGC %d, sf = %2.2f c/deg', iRGC, spatialFrequency));
                xlabel('time (sec)');
                ylabel('response');
                drawnow;
                if (exportFig)
                    pdfFileName = sprintf('Response_RGC_%d_SF_%2.2fcpd_LMS_%0.2f_%0.2f_%0.2f_PolansSID_%d_%s',iRGC, spatialFrequency, LMScontrast(1), LMScontrast(2), LMScontrast(3), PolansSubjectID, opticsPostFix);
                    plotlabOBJ.exportFig(hFig, 'pdf', pdfFileName, exportDir);
                end
                
            end
        end
        
        if (visualizeIndividualFits)
            setupPlotLab(-1);
        end
    
end

function plotlabOBJ = setupPlotLab(mode)
    if (mode == 0)
        plotlabOBJ = plotlab();
        plotlabOBJ.applyRecipe(...
                'colorOrder', [1 0 0; 0 0 1], ...
                'axesBox', 'off', ...
                'axesTickDir', 'both', ...
                'renderer', 'painters', ...
                'lineMarkerSize', 14, ...
                'axesTickLength', [0.01 0.01], ...
                'legendLocation', 'SouthWest', ...
                'figureWidthInches', 7, ...
                'figureHeightInches', 7);
    else
        plotlab.resetAllDefaults();
    end
end 
