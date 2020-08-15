function performPSFConvoComputations(varargin)

% Polans et al subjects grouped according to different criteria
%     sharpestPSFSubjectIDs = [4 9];
%     mediumSharpnessPSFSubjectIDs = [5 8 10];
%     blurriestPSFSubjectIDs = [7];
%     noArtifactPSFSubjectIDs = [4 5 7 8 9 10];
%     someArtifactPSFSubjectIDs = [1 3 6];
%     largeArtifacPSFSubjectIDs = [2];


    % Parse input
    p = inputParser;
    p.addParameter('PolansSubjectIDs', [4], @isnumeric);
    p.addParameter('eccTested', [0 0.2 0.5 1 1.5 2 2.5 3:25]);
    p.addParameter('quadrantsToCompute', {'horizontal'}); %, @(x)(ismember(x, {'horizontal', 'superior', 'inferior'})));
    p.addParameter('generateNewDeconvolutionFiles', false, @islogical);
    p.addParameter('visualizeFits', false, @islogical);
    p.parse(varargin{:});
    
    PolansSubjectIDs = p.Results.PolansSubjectIDs;
    eccTested = p.Results.eccTested;
    quadrantsToCompute = p.Results.quadrantsToCompute;
    generateNewDeconvolutionFiles = p.Results.generateNewDeconvolutionFiles;
    visualizeFits = p.Results.visualizeFits;
    
    
    deconvolutionOpticsParams.PolansWavefrontAberrationSubjectIDsToCompute = PolansSubjectIDs;
    deconvolutionOpticsParams.quadrantsToCompute = quadrantsToCompute;
        
    % Perform the deconvolution analysis for certain Polans subjects 
    if (generateNewDeconvolutionFiles)
      
        % Instantiate a CronerKaplanRGCModel
        ck = CronerKaplanRGCModel(...
            'generateAllFigures', false, ...
            'instantiatePlotLab', false);
        
        % Tell it to generate deconvolution files for the desired
        % eccentricities, and deconvolution parameters (quadrants&subject)
        ck.generateDeconvolutionFiles(deconvolutionOpticsParams, ...
            'eccTested', eccTested, ...        
            'visualizeFits', visualizeFits ...
            );
    end
    
    
    performTests = true;
    if (performTests)
        % Test1: Visualize the deconvolution model
        visualizeDeconvolutionModel(deconvolutionOpticsParams, eccTested);
        
        % Test2: Synthesize RF params
        %synthesizeRFparams(deconvolutionOpticsParams);
    end
end

function visualizeDeconvolutionModel(deconvolutionOpticsParams, eccTested)

    % Instantiate a CronerKaplanRGCModel
    ck = CronerKaplanRGCModel(...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
    
    % Assemble and plot the deconvolution model
    deconvolutionModel = ck.computeDeconvolutionModel(deconvolutionOpticsParams,  eccTested);
    
    ck.plotDeconvolutionModel(deconvolutionModel);
end

function synthesizeRFparams(deconvolutionOpticsParams)

     % assume 1 cone input across all eccentricities
    rgcsNum = 30;
    rfCenterInputConesNum = ones(1, rgcsNum); 
    
    minEccDegs = 0.1;
    maxEccDegs = 30;
    minEccMicrons = WatsonRGCModel.rhoDegsToMMs(minEccDegs)*1e3;
    maxEccMicrons = WatsonRGCModel.rhoDegsToMMs(maxEccDegs)*1e3;
    
    % RF positions
    rfCenterPositionMicrons(:,1) = logspace(log10(minEccMicrons), log10(maxEccMicrons), rgcsNum);
    rfCenterPositionMicrons(:,2) = rfCenterPositionMicrons(:,1) * 0;
    
    % Synthesize params
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(...
        rfCenterInputConesNum, rfCenterPositionMicrons, deconvolutionOpticsParams);
    
    % Plot synthesized params
    figure(3);
    subplot(2,2,1);
    plot(synthesizedRFParams.rfEccRadiusDegs,  synthesizedRFParams.visual.centerCharacteristicRadiiDegs, ...
        'o-', 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]); hold on;
    plot(synthesizedRFParams.rfEccRadiusDegs,  synthesizedRFParams.visual.surroundCharacteristicRadiiDegs, ...
        'o-', 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1]);
    set(gca, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100]);
    set(gca, 'YScale', 'log', 'YLim', [0.003 10], 'YTick', [0.003 0.01 0.03 0.1 0.3 1 3 10]);
    axis 'square';
    xlabel('ecc (degs)');
    ylabel('radius (degs)');
    
    subplot(2,2,2);
    yyaxis 'left'
    plot(synthesizedRFParams.visual.centerCharacteristicRadiiDegs, synthesizedRFParams.visual.centerPeakSensitivities, ...
        'o-', 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]); hold on;
    plot(synthesizedRFParams.visual.centerCharacteristicRadiiDegs, synthesizedRFParams.visual.surroundPeakSensitivities, ...
        'o-', 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1]);
    set(gca, 'YScale', 'log', 'YLim', [0.003 3000], 'YTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000]);
    ylabel('peak sensitivity');
    
    yyaxis 'right'
    plot(synthesizedRFParams.visual.centerCharacteristicRadiiDegs, ...
        synthesizedRFParams.visual.surroundPeakSensitivities ./ synthesizedRFParams.visual.centerPeakSensitivities, '--');
    set(gca, 'YScale', 'log', 'YLim', [0.001 1], 'YTick', [0.001 0.003 0.01 0.03 0.1 0.3 1]);
    set(gca, 'XScale', 'log', 'XLim', [0.003 3], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3]);
    ylabel('surround:center peak sensitivity');
    axis 'square';
    xlabel('radius (degs)');
    
    
    subplot(2,2,3);
    plot(synthesizedRFParams.rfEccRadiusDegs, synthesizedRFParams.visual.centerCharacteristicRadiiDegs ./ synthesizedRFParams.visual.surroundCharacteristicRadiiDegs, 'ko-');
    set(gca, 'XLim', [0 100], 'XTick', 0:10:100);
    set(gca,  'YLim', [0 1], 'YTick', 0:0.1:1);
    axis 'square';
    xlabel('ecc (degs)');
    ylabel('center/surround radius ratio');
    
    subplot(2,2,4);
    scIntSensRatio = (synthesizedRFParams.visual.surroundPeakSensitivities./synthesizedRFParams.visual.centerPeakSensitivities) .* ...
                     (synthesizedRFParams.visual.surroundCharacteristicRadiiDegs./synthesizedRFParams.visual.centerCharacteristicRadiiDegs).^2;
    plot(synthesizedRFParams.rfEccRadiusDegs, scIntSensRatio, 'ko-');
    set(gca, 'XLim', [0 100], 'XTick', 0:10:100);
    set(gca, 'YLim', [0 1], 'YTick', 0:0.1:1);
    xlabel('ecc (degs)');
    ylabel('surround/center integrated sensitivity ratio');
end
