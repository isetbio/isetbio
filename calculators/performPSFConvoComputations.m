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
    p.addParameter('eccTested', [0 0.25 0.5 1 1.5 2.0 2.5 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25]);
    p.addParameter('quadrantsToCompute', {'horizontal'}); %, @(x)(ismember(x, {'horizontal', 'superior', 'inferior'})));
    p.addParameter('generateNewDeconvolutionFiles', false, @islogical);
    p.parse(varargin{:});
    
    PolansSubjectIDs = p.Results.PolansSubjectIDs;
    eccTested = p.Results.eccTested;
    quadrantsToCompute = p.Results.quadrantsToCompute;
    generateNewDeconvolutionFiles = p.Results.generateNewDeconvolutionFiles;
    
    ck = CronerKaplanRGCModel(...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
     
    % Perform the deconvolution analysis for certain Polans subjects 
    if (generateNewDeconvolutionFiles)
        deconvolutionOpticsParams = struct(...
            'PolansWavefrontAberrationSubjectIDsToCompute', PolansSubjectIDs ...
            );
        deconvolutionOpticsParams.quadrantsToCompute = quadrantsToCompute;
        
        ck.generateDeconvolutionFiles(...
            deconvolutionOpticsParams, ...
            'eccTested', eccTested);
    end
    
    
    performTests = ~true;
    if (performTests)
        performDeconvolutionTests(PolansSubjectIDs, quadrantsToCompute)
    end
end

function performDeconvolutionTests(PolansSubjectIDs, quadrantsToCompute)

    plotDeconvolutionModel =  true;
    ck = CronerKaplanRGCModel(...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
    
    % TEST1. Compute and plot deconvolution model for only the 'horizontal'meridian
    deconvolutionOpticsParams = struct(...
        'PolansWavefrontAberrationSubjectIDsToAverage', PolansSubjectIDs ...
    );
    deconvolutionOpticsParams.quadrantsToAverage = quadrantsToCompute;
  
    deconvolutionModel = ck.computeDeconvolutionModel(deconvolutionOpticsParams);
    
    if (plotDeconvolutionModel)
        CronerKaplanRGCModel.plotDeconvolutionModel(deconvolutionModel);
    end
    
    % TEST2. Synthesize params 
    rgcsNum = 30;
    rfCenterInputConesNum = ones(1, rgcsNum);  % assume 1 cone input across all eccentricities
    
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
    plot(synthesizedRFParams.eccDegs,  synthesizedRFParams.visual.centerCharacteristicRadiiDegs, ...
        'o-', 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]); hold on;
    plot(synthesizedRFParams.eccDegs,  synthesizedRFParams.visual.surroundCharacteristicRadiiDegs, ...
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
    plot(synthesizedRFParams.eccDegs, synthesizedRFParams.visual.centerCharacteristicRadiiDegs ./ synthesizedRFParams.visual.surroundCharacteristicRadiiDegs, 'ko-');
    set(gca, 'XLim', [0 100], 'XTick', 0:10:100);
    set(gca,  'YLim', [0 1], 'YTick', 0:0.1:1);
    axis 'square';
    xlabel('ecc (degs)');
    ylabel('center/surround radius ratio');
    
    subplot(2,2,4);
    scIntSensRatio = (synthesizedRFParams.visual.surroundPeakSensitivities./synthesizedRFParams.visual.centerPeakSensitivities) .* ...
                     (synthesizedRFParams.visual.surroundCharacteristicRadiiDegs./synthesizedRFParams.visual.centerCharacteristicRadiiDegs).^2;
    plot(synthesizedRFParams.eccDegs, scIntSensRatio, 'ko-');
    set(gca, 'XLim', [0 100], 'XTick', 0:10:100);
    set(gca, 'YLim', [0 1], 'YTick', 0:0.1:1);
    xlabel('ecc (degs)');
    ylabel('surround/center integrated sensitivity ratio');
end
