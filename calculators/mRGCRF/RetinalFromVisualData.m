function RetinalFromVisualData()
    % Step 1. Convolve different size Gaussians with good subject (Polans et al) PSFs
    % convolveGaussianWithPSF()
    
    % Step 2. Analyze the results of this convolution
    analyzeConvolutionResults()
    
    % Step 3. Plot the corrected visual data based on the convolution  results
    plotCorrectedCronerKaplanStats()
end

