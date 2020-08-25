function generateDeconvolutionFilesForMidgetRGCs(varargin)
% Generate the deconvolution files for a particular Polans subject and
% eccentricity quadrant
%
% Syntax:
%   generateDeconvolutionFilesForMidgetRGCs(varargin)
%
% Description:
%    Generates the various deconvolution files which are used to compute the 
%    characteristic radii and peak sensitivities for retinal/visual mRGC RFs.
%
% Outputs:
%    None. Deconvolution files are written in  isetbio/calculators/mRGCRF/VisualToRetinalCorrectionData/DecolvolutionData
%
% Optional key/value pairs:
%    PolansSubjectIDs               : [1 - 10]. See comments below about the different subjects
%    eccTested                      : Vector of eccentricities at which to compute the deconvolution data
%    examinedConesNumInRFCenter     : Vector of number of cones in the RF center
%    quadrantsToCompute             : 'horizontal', 'superior', 'inferior'
%    generateNewDeconvolutionFiles  : whether to generate the files or just  visualize stuff
%    visualizeFits                  : whether to visualize the various fits
%
% Notes about the Polans et al subjects:
%     sharpestPSFSubjectIDs = [4 8 9];
%     mediumSharpnessPSFSubjectIDs = [5 10];
%     blurriestPSFSubjectIDs = [7];
%     noArtifactPSFSubjectIDs = [4 5 7 8 9 10];
%     someArtifactPSFSubjectIDs = [1 3 6];
%     largeArtifacPSFSubjectIDs = [2];
% 
% Example usages:
%{
   % Generate figures for presentation:

   generateDeconvolutionFilesForMidgetRGCs('PolansSubjectID', 4, ...
        'subregion', 'center', ...
        'eccTested', 1.5, ...
        'examinedConesNumInRFCenter', 1, 'generateNewDeconvolutionFiles', true, ...
        'visualizeFits', true, 'exportFig', true, ...
        'visualizeDeconvolutionModel', false, ...
        'synthesizeDemoRFparams', false);


   % Generate CENTER deconvolution data for subject 4 optics and the default ecc
     range (i.e., [0 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2 2.5 3:20]) and for
     RF centers containing between 1 and 20 cones, skipping visualization.

   generateDeconvolutionFilesForMidgetRGCs('PolansSubjectID', 4, ...
        'subregion', 'center', ...
        'examinedConesNumInRFCenter', 1:20, 'generateNewDeconvolutionFiles', true, ...
        'visualizeFits', ~true, 'exportFig', ~true, ...
        'visualizeDeconvolutionModel', false, ...
        'synthesizeDemoRFparams', false);


    % Generate SURROUND deconvolution data. It depends on the center
    deconvolution data so we run after the CENTER deconvolution data are
    available. Also it takes much longer to compute, so we run on two
    separate MATLAB instances by splitting the eccentricities in even/odd:

    % All eccentricities:
    eccTested =  [0 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2 2.5 3:20];

    INSTANCE-1: one half of eccentricities
    generateDeconvolutionFilesForMidgetRGCs('PolansSubjectID', 4, ...
        'subregion', 'surround', ...
        'eccTested', eccTested(1:2:end), ...
        'examinedConesNumInRFCenter', 1:20, ...
        'generateNewDeconvolutionFiles', true, ...
        'visualizeFits', ~true, 'exportFig', ~true, ...
        'visualizeDeconvolutionModel', false, ...
        'synthesizeDemoRFparams', false);

    INSTANCE-2: second half of eccentricities
    generateDeconvolutionFilesForMidgetRGCs('PolansSubjectID', 4, ...
        'subregion', 'surround', ...
        'eccTested', eccTested(1:2:end), ...
        'examinedConesNumInRFCenter', 1:20, ...
        'generateNewDeconvolutionFiles', true, ...
        'visualizeFits', ~true, 'exportFig', ~true, ...
        'visualizeDeconvolutionModel', false, ...
        'synthesizeDemoRFparams', false);

%}

% History:
%    8/2020  NPC  ISETBIO TEAM, 2020

    % Parse the input
    p = inputParser;
    p.addParameter('PolansSubjectIDs', [4], @isnumeric);
    p.addParameter('eccTested', [0 0.1 0.2 0.3 0.4 0.5 0.7 1 1.5 2 2.5 3:20] );
    p.addParameter('examinedConesNumInRFCenter', 1:30);
    p.addParameter('subregion', 'surround', @(x)(ismember(x, {'center', 'surround'})));
    p.addParameter('quadrantsToCompute', {'horizontal'}); %, @(x)(ismember(x, {'horizontal', 'superior', 'inferior'})));
    p.addParameter('generateNewDeconvolutionFiles', false, @islogical);
    p.addParameter('visualizeFits', false, @islogical);
    p.addParameter('exportFig', false, @islogical);
    p.addParameter('visualizeDeconvolutionModel', false, @islogical);
    p.addParameter('synthesizeDemoRFparams', true, @islogical);
    p.parse(varargin{:});
    
    PolansSubjectIDs = p.Results.PolansSubjectIDs;
    eccTested = p.Results.eccTested;
    examinedConesNumInRFCenter = p.Results.examinedConesNumInRFCenter;
    subregion = p.Results.subregion;
    quadrantsToCompute = p.Results.quadrantsToCompute;
    generateNewDeconvolutionFiles = p.Results.generateNewDeconvolutionFiles;
    visualizeFits = p.Results.visualizeFits;
    exportFig = p.Results.exportFig;
    synthesizeDemoRFparams = p.Results.synthesizeDemoRFparams;
    visualizeDeconvolutionModel = p.Results.visualizeDeconvolutionModel;
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
        ck.generateDeconvolutionFiles(deconvolutionOpticsParams, subregion, ...  
            'examinedConesNumInRFCenter', examinedConesNumInRFCenter, ...
            'visualizeFits', visualizeFits, ...
            'exportFig', exportFig...
            );
    end
    
    if (visualizeDeconvolutionModel)
        % Compute and visualize the deconvolution model 
        computeAndVisualizeDeconvolutionModel(deconvolutionOpticsParams);
    end
    
    if (synthesizeDemoRFparams)
        % Synthesize RF params for single cone input RF centers at different eccentricities
        synthesizeRFparams(deconvolutionOpticsParams);
    end
end

function computeAndVisualizeDeconvolutionModel(deconvolutionOpticsParams)
    % Instantiate a CronerKaplanRGCModel
    ck = CronerKaplanRGCModel(...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
    
    % Assemble and plot the deconvolution model
    deconvolutionModel = ck.computeDeconvolutionModel(deconvolutionOpticsParams);
    
    ck.plotDeconvolutionModel(deconvolutionModel);
end

function synthesizeRFparams(deconvolutionOpticsParams)

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
    maxEccDegs = 5.0;
    minEccMicrons = WatsonRGCModel.rhoDegsToMMs(minEccDegs)*1e3;
    maxEccMicrons = WatsonRGCModel.rhoDegsToMMs(maxEccDegs)*1e3;
    
    % RF positions
    rfCenterPositionMicrons(:,1) = logspace(log10(minEccMicrons), log10(maxEccMicrons), rgcsNum);
    rfCenterPositionMicrons(:,2) = rfCenterPositionMicrons(:,1) * 0;
    
    % Synthesize params
    synthesizedRFParams = ck.synthesizeRetinalRFparamsConsistentWithVisualRFparams(...
        rfCenterInputConesNum, rfCenterPositionMicrons, deconvolutionOpticsParams);
    
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
