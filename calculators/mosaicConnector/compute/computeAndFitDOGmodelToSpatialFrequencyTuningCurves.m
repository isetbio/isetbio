function [patchDogParams, responseAmplitude, spatialFrequenciesCPDHR, responseAmplitudeHR, responseTimeAxisHR, fittedResponsesHR] = ...
    computeAndFitDOGmodelToSpatialFrequencyTuningCurves(...
    responseTimeAxis, integratedResponsesMean, integratedResponsesStDev, ...
    maxSpikeRateModulation, stimSpatialParams, stimTemporalParams,  spatialFrequenciesCPD, visualizeIndividualFits, ...
    targetRGCs, LMScontrast,  opticsPostFix, PolansSubjectID, exportFig, figExportsDir, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('synthParams', [], @(x)(isempty(x)||(isstruct(x))));
    p.parse(varargin{:});
    synthParams = p.Results.synthParams;
    
    visualizeResponseTimeCourseFits = false;
    
    % Compute response modulations at each spatial frequency
    for sfIndex = 1:numel(spatialFrequenciesCPD)
        switch (stimSpatialParams.type)
            case 'driftingGrating'
                visualizeIndividualSinusoidalFits = visualizeResponseTimeCourseFits;
                exportSinusoidalFig = visualizeResponseTimeCourseFits;
                % Compute modulation of the response at the fundamental temporal frequency
                [responseAmplitude(:, sfIndex), ...
                 responseAmplitudeSE(:, sfIndex), ...
                 responsePhase(:, sfIndex), ...
                 responseTimeAxisHR, fittedResponsesHR(sfIndex,:,:)] = fitSinusoidalModulationToResponseTimeCourse(...
                    squeeze(integratedResponsesMean(sfIndex,:,:)), ...
                    squeeze(integratedResponsesStDev(sfIndex,:,:)), ...
                    responseTimeAxis, ...
                    stimTemporalParams.temporalFrequencyHz, ...
                    spatialFrequenciesCPD(sfIndex), maxSpikeRateModulation, ...
                     visualizeIndividualSinusoidalFits, targetRGCs, LMScontrast,  opticsPostFix, PolansSubjectID, exportSinusoidalFig, figExportsDir);
            otherwise
                error('Unknown stimulus type: ''%''.', stimulusSpatialParams.type)
        end
        
    end % sfIndex
    
    % Fit the DoG model to the SF tuning curves
    DoGFunction = @(params,sf)(...
            params(1)           * ( pi * (params(2))^2 * exp(-(pi*params(2)*sf).^2) ) - ...
            params(1)*params(3) * ( pi * (params(2)*params(4))^2 * exp(-(pi*params(2)*params(4)*sf).^2) ) ...
    ); 

    %                 Kc   Rc     kS/kC  Rs/Rc
    initialParams = [300   0.05    1e-2   7];
    
    % Upper and lower values of DoG params
    %               Kc       Rc     kS/kC       Rs/Rc  
    lowerBounds   = [1     0.0001    1e-5         2  ];
    upperBounds   = [Inf   2.0       1e0          20 ];
    
    [patchDogParams,spatialFrequenciesCPDHR, responseAmplitudeHR] = ...
        fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseAmplitude, responseAmplitudeSE, ...
        DoGFunction, initialParams, lowerBounds, upperBounds, ...
        maxSpikeRateModulation, visualizeIndividualFits, ...
        LMScontrast, opticsPostFix, PolansSubjectID, exportFig, figExportsDir, ...
        'synthParams', synthParams);
end


function [patchDogParams, spatialFrequenciesCPDHR, responseTuningHR] = ...
    fitDoGmodelToSpatialFrequencyCurve(spatialFrequenciesCPD, responseTuning, responseTuningSE,...
    DoGFunction, initialParams,  lowerBounds, upperBounds, ......
    maxSpikeRateModulation, visualizeIndividualFits, ...
    LMScontrast, opticsPostFix, PolansSubjectID, exportFig, exportDir, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('synthParams', [], @(x)(isempty(x)||(isstruct(x))));
    p.parse(varargin{:});
    synthParams = p.Results.synthParams;
    

    % Fitting options
    oldoptions = optimoptions('lsqcurvefit');
    options = optimoptions(oldoptions,'MaxFunctionEvaluations',5000, 'FunctionTolerance', 1e-8, 'MaxIterations', 10000);
    
    % Preallocate memory
    rgcsNum = size(responseTuning,1);
    patchDogParams = cell(1, rgcsNum);
    fittedParamsAllRGCs = zeros(rgcsNum, length(lowerBounds));
    
    % High-resolution spatial frequency axis
    spatialFrequenciesCPDHR = logspace(log10(spatialFrequenciesCPD(1)), log10(spatialFrequenciesCPD(end)), 100);
    responseTuningHR = zeros(rgcsNum, numel(spatialFrequenciesCPDHR));
    
    % Fit spatial frequency curves for each RGC
    for iRGC = 1:rgcsNum
        % Retrieve data to be fitted
        theSFtuning = responseTuning(iRGC,:);
        theSFtuningSE = responseTuningSE(iRGC,:);
        
        % Fit the model to the data
        fittedParams = lsqcurvefit(DoGFunction, initialParams, spatialFrequenciesCPD, theSFtuning,lowerBounds,upperBounds, options);
        
        % Do a multi-start to find global minimum
        problem = createOptimProblem('lsqcurvefit',...
            'x0',fittedParams, ...
            'objective',DoGFunction,...
            'lb',lowerBounds, ...
            'ub',upperBounds,...
            'xdata',spatialFrequenciesCPD,...
            'ydata',theSFtuning);

        displayProgress = 'off'; % 'iter';
        ms = MultiStart(...
            'Display', displayProgress, ...
            'FunctionTolerance', 2e-4, ...
            'UseParallel', true);

        [fittedParams,errormulti] = run(ms,problem,50);

        patchDogParams{iRGC} = struct(...
            'kC', fittedParams(1), ...
            'rC', fittedParams(2), ...
            'kS', fittedParams(1)*fittedParams(3), ...
            'rS', fittedParams(2)*fittedParams(4) ...
            );
        
        
        % Compute model fit on a high-res spatial frequency axis
        responseTuningHR(iRGC,:) = DoGFunction(fittedParams, spatialFrequenciesCPDHR);
        
        % Visualize model fit and data
        if (visualizeIndividualFits)
            visualizeSpatialFrequencyTuning([], spatialFrequenciesCPD, theSFtuning, theSFtuningSE, maxSpikeRateModulation, ...
                spatialFrequenciesCPDHR, squeeze(responseTuningHR(iRGC,:)), patchDogParams{iRGC}, ...
                iRGC, LMScontrast, opticsPostFix, PolansSubjectID, exportFig, exportDir, ...
                'synthParams', synthParams);
        end

    end % iRGC
    
end

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
            
            if ((visualizeIndividualFits) && (ismember(iRGC, targetRGC)))
                hFig = figure(333); clf;
                errorbar(responseTimeAxis, response, 1./weights, 'ko-', 'LineWidth', 1.0, 'Color', [0.4 0.4 0.4]); hold on;
                line(responseTimeAxisHR, responseHR, 'Color', [0 1 1], 'LineWidth', 3); hold off;
                set(gca, 'YLim', maxSpikeRate*[-1 1], 'XLim', [responseTimeAxis(1) responseTimeAxis(end)]);
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
