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
        deconvolutionOpticsParams.quadrantsToCompute =  quadrantsToCompute;
        
        ck.generateDeconvolutionFiles(...
            deconvolutionOpticsParams, ...
            'eccTested', eccTested);
    end
    
    
    % Compute and plot deconvolution model for only the 'horizontal'meridian
    deconvolutionOpticsParams = struct(...
        'PolansWavefrontAberrationSubjectIDsToAverage', PolansSubjectIDs ...
    );
    deconvolutionOpticsParams.quadrantsToAverage = {'horizontal'};
  
    deconvolutionModel = ck.computeDeconvolutionModel(deconvolutionOpticsParams);
    CronerKaplanRGCModel.plotDeconvolutionModel(deconvolutionModel);
    
end
