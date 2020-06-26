function performPSFConvoComputations

    ck = CronerKaplanRGCModel(...
        'generateAllFigures', false, ...
        'instantiatePlotLab', false);
    
    ck.performGaussianConvolutionWithPolansPSFanalysis(...
        'eccTested', 20:25);

end
