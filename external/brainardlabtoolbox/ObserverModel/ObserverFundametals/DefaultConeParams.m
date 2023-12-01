function coneParams = DefaultConeParams(type)
% Generate default cone params structure
%
% Syntax:
%    coneParams = DefaultConeParams(type)
%
% Description:
%    Generate a structure describing and observer's cone fundamentals, with
%    reasonable defaults.
%   
%    The input string type allows some flexibility about the description.
%
%    Asano et al. give the following population SD's for the individual
%    difference parameters (their Table 5, Step 2 numbers:
%       Lens    - 18.7%
%       Macular - 36.5%
%       L Density - 9%
%       M Density - 9%
%       S Density - 7.4%
%       L Shift   - 2 nm
%       M Shift   - 1.5 nm
%       S Shift   - 1.3 nm
%
% Inputs:
%     type                          - String specifying cone parameterization type.
%                                     'cie_asano': The CIE fundamentals
%                                      with Asano et al. individual
%                                      differene paramters.
%
% Outputs:
%     coneParams                    - Structure with field for each parameter.
%
% Optional key/value pairs:
%     None.
%
% See also: ObserverVecToParams, ObserverParamsToVec
%

% History:
%   08/10/19  dhb  Wrote it.

% Examples:
%{
    coneParams = DefaultConeParams('cie_asano')
    coneParams.indDiffParams
%}

switch (type)
    case 'cie_asano'
        
        coneParams.type = 'cie_asano';
        
        % Basic CIE parameters
        coneParams.fieldSizeDegrees = 10;
        coneParams.ageYears = 32;
        coneParams.pupilDiamMM = 3;
        
        % Asano individual difference params
        coneParams.indDiffParams.dlens = 0;
        coneParams.indDiffParams.dmac = 0;
        coneParams.indDiffParams.dphotopigment = [0 0 0]';
        coneParams.indDiffParams.lambdaMaxShift = [0 0 0]';
        coneParams.indDiffParams.shiftType = 'linear';
        
    otherwise
        error('Unknown cone parameters type passed.');
end

