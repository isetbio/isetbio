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

