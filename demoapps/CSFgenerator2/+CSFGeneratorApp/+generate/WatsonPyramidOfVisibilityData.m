function [sfSupport, S] = WatsonPyramidOfVisibilityData(sfMin, sfMax, meanLuminanceCdM2, constantParam)
    temporalFrequency = 0.0; cW = 0;
    logLuminanceNits = log10(meanLuminanceCdM2);
    % Table 1, of Watson 2018, "The Field of View, the Field of Resolution and the 
    %       Field of Contrast Sensitivity" (Luminance)
    if (strcmp(constantParam, 'constant cycles'))
        c0 = 1.380; cF = -0.091; cL = 0.391;
    else
        c0 = 1.739; cF = -0.060; cL = 0.391;
    end
    
    sfSupport = sfMin:1:sfMax;
    logS = c0 + cW*temporalFrequency + cF*sfSupport + cL*logLuminanceNits;
    S = 10.^logS;
end