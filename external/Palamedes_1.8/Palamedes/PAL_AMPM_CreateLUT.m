%
%PAL_AMPM_CreateLUT  Look-up Table (LUT) of Psychometric Function values
%
%   syntax: PFLookUpTable = PAL_AMPM_CreateLUT(priorAlphaValues, ...
%       priorBetaValues, priorGammaValues, priorLambdaValues, ...
%       StimLevels, PF, gammaEQlambda)
%
%   Creates a 5-D (alpha x beta x gamma x lambda x stimulus intensity) 
%   array of PF values.
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.2.0, 1.5.0, 1.6.0, 1.6.3 (see History.m)

function PFLookUpTable = PAL_AMPM_CreateLUT(priorAlphaValues, priorBetaValues, priorGammaValues, priorLambdaValues, StimLevels, PF, gammaEQlambda)
    
[a, b, g, l, x] = ndgrid(priorAlphaValues, priorBetaValues, priorGammaValues, priorLambdaValues, StimLevels);
params.alpha = a;
params.beta = 10.^b;
params.gamma = g;
params.lambda = l;

if gammaEQlambda
    params.gamma = l;
end

PFLookUpTable = PF(params, x);