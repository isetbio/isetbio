function generateDeconvolutionFilesForMidgetRGCs(varargin)

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
            'deconvolutionEccs', eccTested, ...
            'generateAllFigures', false, ...
            'instantiatePlotLab', false);
        
        % Tell it to generate deconvolution files for the desired
        % eccentricities, and deconvolution parameters (quadrants&subject)
        ck.generateDeconvolutionFiles(deconvolutionOpticsParams, ...      
            'visualizeFits', visualizeFits ...
            );
    end
    
    
    visualizeDeconvolutionModel(deconvolutionOpticsParams, eccTested);
    
    performTests = true;
    if (performTests)
        % Test2: Synthesize RF params
        synthesizeRFparams(deconvolutionOpticsParams, eccTested);
    end
end

function visualizeDeconvolutionModel(deconvolutionOpticsParams, eccTested)

    % Instantiate a CronerKaplanRGCModel
    ck = CronerKaplanRGCModel(...
        'deconvolutionEccs', eccTested, ...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
    
    % Assemble and plot the deconvolution model
    deconvolutionModel = ck.computeDeconvolutionModel(deconvolutionOpticsParams);
    
    ck.plotDeconvolutionModel(deconvolutionModel);
end

function synthesizeRFparams(deconvolutionOpticsParams, eccTested)

    % Instantiate a CronerKaplanRGCModel
    ck = CronerKaplanRGCModel(...
            'generateAllFigures', false, ...
            'instantiatePlotLab', false);
        
    % Number of RGCs to synthesize
    rgcsNum = 30;
    
    % assume 1 cone input across all eccentricities
    conesInRFCenter = 1;
    rfCenterInputConesNum = ones(1, rgcsNum)*conesInRFCenter; 
    
    minEccDegs = 0.01;
    maxEccDegs = 3.0;
    minEccMicrons = WatsonRGCModel.rhoDegsToMMs(minEccDegs)*1e3;
    maxEccMicrons = WatsonRGCModel.rhoDegsToMMs(maxEccDegs)*1e3;
    
    % RF positions
    rfCenterPositionMicrons(:,1) = logspace(log10(minEccMicrons), log10(maxEccMicrons), rgcsNum);
    rfCenterPositionMicrons(:,2) = rfCenterPositionMicrons(:,1) * 0;
    
    % Synthesize params
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(...
        rfCenterInputConesNum, rfCenterPositionMicrons, deconvolutionOpticsParams, eccTested);
    
    % Plot synthesized params
    figNo = 3;
    plotSynthesizedParams(figNo, synthesizedRFParams.rfEccRadiusDegs, synthesizedRFParams.visual, 'visual');
    
    figNo = 4;
    plotSynthesizedParams(figNo, synthesizedRFParams.rfEccRadiusDegs, synthesizedRFParams.retinal, 'retinal');
    
end

function plotSynthesizedParams(figNo, rfEccRadiusDegs, synthesizedRFParams, domain)
    figure(figNo); clf;
    subplot(2,2,1);
    plot(rfEccRadiusDegs,  synthesizedRFParams.centerCharacteristicRadiiDegs, ...
        'o-', 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]); hold on;
    plot(rfEccRadiusDegs,  synthesizedRFParams.surroundCharacteristicRadiiDegs, ...
        'o-', 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1]);
    set(gca, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100]);
    set(gca, 'YScale', 'log', 'YLim', [0.003 10], 'YTick', [0.003 0.01 0.03 0.1 0.3 1 3 10]);
    axis 'square';
    xlabel('ecc (degs)');
    ylabel(sprintf('%s radius (degs)', domain));
    
    subplot(2,2,2);
    yyaxis 'left'
    plot(synthesizedRFParams.centerCharacteristicRadiiDegs, synthesizedRFParams.centerPeakSensitivities, ...
        'o-', 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [1 0.5 0.5]); hold on;
    plot(synthesizedRFParams.centerCharacteristicRadiiDegs, synthesizedRFParams.surroundPeakSensitivities, ...
        'o-', 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0.5 0.5 1]);
    set(gca, 'YScale', 'log', 'YLim', [0.003 3000], 'YTick', [0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300 1000 3000]);
    ylabel('peak sensitivity');
    
    yyaxis 'right'
    plot(synthesizedRFParams.centerCharacteristicRadiiDegs, ...
        synthesizedRFParams.surroundPeakSensitivities ./ synthesizedRFParams.centerPeakSensitivities, '--');
    set(gca, 'YScale', 'log', 'YLim', [0.001 1], 'YTick', [0.001 0.003 0.01 0.03 0.1 0.3 1]);
    set(gca, 'XScale', 'log', 'XLim', [0.003 3], 'XTick', [0.003 0.01 0.03 0.1 0.3 1 3]);
    axis 'square';
    ylabel(sprintf('%s surround:center peak sensitivity', domain));
    
    
    subplot(2,2,3);
    plot(rfEccRadiusDegs, synthesizedRFParams.centerCharacteristicRadiiDegs ./ synthesizedRFParams.surroundCharacteristicRadiiDegs, 'ko-');
    set(gca, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100]);
    set(gca,  'YLim', [0 1], 'YTick', 0:0.1:1);
    xlabel('ecc (degs)');
    ylabel(sprintf('%s center/surround radius ratio', domain));
    
    subplot(2,2,4);
    scIntSensRatio = (synthesizedRFParams.surroundPeakSensitivities./synthesizedRFParams.centerPeakSensitivities) .* ...
                     (synthesizedRFParams.surroundCharacteristicRadiiDegs./synthesizedRFParams.centerCharacteristicRadiiDegs).^2;
    plot(rfEccRadiusDegs, scIntSensRatio, 'ko-');
    set(gca, 'XScale', 'log', 'XLim', [0.1 100], 'XTick', [0.1 0.3 1 3 10 30 100]);
    set(gca, 'YLim', [0 1], 'YTick', 0:0.1:1);
    xlabel('ecc (degs)');
    ylabel(sprintf('%s surround/center integrated sensitivity ratio', domain));
end
