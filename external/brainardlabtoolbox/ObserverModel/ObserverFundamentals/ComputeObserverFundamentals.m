function [T,T_energy,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ComputeObserverFundamentals(coneParams,S)
% Compute cone fundamentals from cone parameter structure
%
% Syntax:
%     T = ComputeObserverFundamentals(coneParams,S)
%
% Description:
%     Compute cone fundamentals from a structure describing the parameters
%     of a cone fundamentlals model.
%
%     T matrix is fundamentals are in energy units and normalized to a max of one.
%     Can also pick up not normalized in either energy or quantal units as
%     additional returns.
%
% Inputs:
%     coneParams                      - Structure providing parameters for
%                                       the fundamentals.  See
%                                       DefaultConeParams.
%     S                               - Wavelength sampling, PTB conventions.
% 
% Outputs:
%     T                               - Matrix of fundamentals, PTB matrix
%                                       format with each fundamental in a
%                                       row.  Normalized in energy units.
%     T_energy                  - Matrix of fundamentals, PTB matrix
%                                       format with each fundamental in a
%                                       row.  Not normalized, in energy units.
%     T_quantal                 - Matrix of fundamentals, PTB matrix
%                                       format with each fundamental in a row.
%                                       Not normalized, and in quantal
%                                       units.
%     adjIndDiffParams          - The Asano adjusted individual difference
%                                       parameters returned by
%                                       CIEComputeConeFundamentals.
%     params                    - Parameters structure returned by
%                                       CIEComputeConeFundamentals.
%     staticParams              - Static parameters structure returned by
%                                       CIEComputeConeFundamentals.
%
% Optional key/value pairs:
%    None.
%
% See also: DefaultConeParams, ObserverParamsToVec, ObserverVecToParams
%

% History:
%   08/10/19  dhb  Wrote it.

% Examples:
%{
    coneParams = DefaultConeParams('cie_asano');
    S = [400 1 301];
    T = ComputeObserverFundamentals(coneParams,S);
	figure; clf; hold on;
    plot(SToWls(S),T(1,:)','r','LineWidth',2);
    plot(SToWls(S),T(2,:)','g','LineWidth',2);
    plot(SToWls(S),T(3,:)','b','LineWidth',2);
    xlabel('Wavelength (nm)');
    ylabel('Fundamental');
%}
%{
    coneParams = DefaultConeParams('cie_asano');
    S = [400 1 301];
    [~,~,~,adjIndDiffParams,params,staticParams] = ComputeObserverFundamentals(coneParams,S);	
%}

switch (coneParams.type)
    case 'cie_asano'
        % Get cone spectral sensitivities
        [~,~,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
            ComputeCIEConeFundamentals(MakeItS(S),coneParams.fieldSizeDegrees,coneParams.ageYears,coneParams.pupilDiamMM, ...
            [],[],[], ...
            [],[],[],coneParams.indDiffParams);
        T_energy = EnergyToQuanta(S,T_quantalIsomerizations')';
        for ii = 1:3
            T(ii,:) = T_energy(ii,:)/max(T_energy(ii,:));
        end

    case {'cie_govardovskii', 'cie_dawis', 'cie_baylor', 'cie_lamb', 'cie_stockmansharpe', 'cie_stockmanrider', 'cie_carrollneitz'}
        % Pop in one of the standard photopigment nomograms, along with the other ind
        % difference parameters.
        useLambdaMax = coneParams.lambdaMax + coneParams.indDiffParams.lambdaMaxShift(:);
        useIndDiffParams = coneParams.indDiffParams;
        useIndDiffParams.lambdaMaxShift = zeros(size(coneParams.indDiffParams.lambdaMaxShift(:)));
        [~,~,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
            ComputeCIEConeFundamentals(MakeItS(S),coneParams.fieldSizeDegrees,coneParams.ageYears,coneParams.pupilDiamMM, ...
            useLambdaMax,coneParams.nomogram,[], ...
            [],[],[],useIndDiffParams);
        T_energy = EnergyToQuanta(S,T_quantalIsomerizations')';
        for ii = 1:3
            T(ii,:) = T_energy(ii,:)/max(T_energy(ii,:));
        end

    otherwise
        error('Unknown cone parameters type passed.');
end


