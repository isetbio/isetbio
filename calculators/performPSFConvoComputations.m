function performPSFConvoComputations

    % Polans et al subjects grouped according to different criteria
    sharpestPSFSubjectIDs = [4 9];  % Subjects with the sharpest PSFs
    mediumSharpnessPSFSubjectIDs = [5 8 10];
    blurriestPSFSubjectIDs = [7];
    noArtifactPSFSubjectIDs = [4 5 7 8 9 10];
    someArtifactPSFSubjectIDs = [1 3 6];
    largeArtifacPSFSubjectIDs = [2];
    
    ck = CronerKaplanRGCModel(...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
%     
    deconvolutionOpticsParams = struct(...
        'PolansWavefrontAberrationSubjectIDsToCompute', 8 ...
        );
    deconvolutionOpticsParams.quadrantsToCompute = {'horizontal', 'upper vertical', 'lower vertical'};
     
    % Generate the various 'ecc_%2.1f_deconvolutions_refractionError_%2.2fD.mat' files, located in the DeconvolutionData directory
    ck.performGaussianConvolutionWithPolansPSFanalysis(...
        deconvolutionOpticsParams, ...
        'eccTested', [0]);
    
    % Generate the deconvolution model only for the horizontal meridian
    deconvolutionOpticsParams = struct(...
        'PolansWavefrontAberrationSubjectIDsToAverage', 8 ...
    );
    deconvolutionOpticsParams.quadrantsToAverage = {'horizontal'};
  
    % Generate and plot the deconvolution model
    ck.generateDeconvolutionModel(...
         deconvolutionOpticsParams, ...
         'modelPrefix', 'horizontalQuadrant_subject8');
end
