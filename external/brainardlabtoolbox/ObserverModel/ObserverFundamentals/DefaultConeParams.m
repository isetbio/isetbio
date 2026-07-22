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
%    Asano et al. give the following population SD's for the individual
%    difference parameters (their Table 5, Step 2 numbers):
%       Lens    - 18.7%
%       Macular - 36.5%
%       L Density - 9%
%       M Density - 9%
%       S Density - 7.4%
%       L Shift   - 2 nm
%       M Shift   - 1.5 nm
%       S Shift   - 1.3 nm
%
% See Asano, Y., Fairchild, M. D., & Blondé, L. (2016). Individual
% colorimetric observer model. PloS one, 11(2), e0145671.
%
% Inputs:
%     type                          - String specifying cone parameterization type.
%                                        'cie_asano': The CIE fundamentals
%                                        with Asano et al. individual
%                                        difference paramters.
%                                     Other types allow swapping in of
%                                     various nomograms for the CIE
%                                     absorbance.
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
        coneParams.nomogram = '';
        coneParams.lambdaMax = [];
        
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

    case 'cie_stockmansharpe'
        coneParams.type = 'cie_stockmansharpe';
        coneParams.nomogram = 'StockmanSharpe';
        coneParams.lambdaMax = [558.9 530.3 420.7]';
        
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
    
    case 'cie_baylor'
        % It's possible 430 would be a more appropriate S cone lambda max.
        coneParams.type = 'cie_baylor';
        coneParams.nomogram = 'Baylor';
        coneParams.lambdaMax = [561 531 420]';
        
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

    case 'cie_govardovskii'
        coneParams.type = 'cie_govardovskii';
        coneParams.nomogram = 'Govardovskii';
        coneParams.lambdaMax = [561 531 420]';

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

    case 'cie_dawis'
        coneParams.type = 'cie_dawis';
        coneParams.nomogram = 'Dawis';
        coneParams.lambdaMax = [561 531 420]';

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

    case 'cie_lamb'
        coneParams.type = 'cie_lamb';
        coneParams.nomogram = 'Lamb';
        coneParams.lambdaMax = [561 531 418.72]';

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

    case 'cie_stockmanrider'
        coneParams.type = 'cie_stockmanrider';
        coneParams.nomogram = 'StockmanRider';

        % Default lambdaMax are the SR template peaks [L; M; S] (nm).
        coneParams.lambdaMax = [551.9 529.8 416.9]';

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

    case 'cie_carrollneitz'
        coneParams.type = 'cie_carrollneitz';
        coneParams.nomogram = 'CarrollNeitz';

        % Default lambdaMax are as in CarrollNeitzNomogram header comments.
        coneParams.lambdaMax = [557.5, 530, 420]';

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

%{
close all;

S = [380 1 401];
wls = SToWls(S);

% Stockman-Rider, baseline
coneParams = DefaultConeParams('cie_stockmanrider');
[T,T_energy,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
    ComputeObserverFundamentals(coneParams,S);
srAborbance = adjIndDiffParams.absorbance;

% Stockman-Sharpe
coneParams = DefaultConeParams('cie_stockmansharpe');
[T,T_energy,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
    ComputeObserverFundamentals(coneParams,S);
ssAborbance = adjIndDiffParams.absorbance;

ssFigure = figure; clf; 
subplot(2,3,1); hold on
plot(wls,srAborbance(1,:)','r','LineWidth',4);
title('Stockman-Sharpe/Stockman-Rider');
subplot(2,3,2); hold on
plot(wls,srAborbance(2,:)','g','LineWidth',4);
title('Stockman-Sharpe/Stockman-Rider');
subplot(2,3,3); hold on
plot(wls,srAborbance(3,:)','b','LineWidth',4);
title('Stockman-Sharpe/Stockman-Rider');
subplot(2,3,1); hold on
plot(wls,ssAborbance(1,:)','k','LineWidth',2);
subplot(2,3,2);
plot(wls,ssAborbance(2,:)','k','LineWidth',2);
subplot(2,3,3);
plot(wls,ssAborbance(3,:)','k','LineWidth',2);

subplot(2,3,4); hold on
plot(wls,log10(srAborbance(1,:)'),'r','LineWidth',4);
title('Stockman-Sharpe/Stockman-Rider');
subplot(2,3,5); hold on
plot(wls,log10(srAborbance(2,:)'),'g','LineWidth',4);
title('Stockman-Sharpe/Stockman-Rider');
subplot(2,3,6); hold on
plot(wls,log10(srAborbance(3,:)'),'b','LineWidth',4);
title('Stockman-Sharpe/Stockman-Rider');
subplot(2,3,4); hold on
plot(wls,log10(ssAborbance(1,:)'),'k','LineWidth',2);
subplot(2,3,5);
plot(wls,log10(ssAborbance(2,:)'),'k','LineWidth',2);
subplot(2,3,6);
plot(wls,log10(ssAborbance(3,:)'),'k','LineWidth',2);

% Carroll-Neitz
coneParams = DefaultConeParams('cie_carrollneitz');
[T,T_energy,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
    ComputeObserverFundamentals(coneParams,S);
cnAborbance = adjIndDiffParams.absorbance;

cnFigure = figure; clf; 
subplot(2,3,1); hold on
plot(wls,srAborbance(1,:)','r','LineWidth',4);
title('Carroll-Neitz/Stockman-Rider');
subplot(2,3,2); hold on
plot(wls,srAborbance(2,:)','g','LineWidth',4);
title('Carroll-Neitz/Stockman-Rider');
subplot(2,3,3); hold on
plot(wls,srAborbance(3,:)','b','LineWidth',4);
title('Carroll-Neitz/Stockman-Rider');
subplot(2,3,1); hold on
plot(wls,cnAborbance(1,:)','k','LineWidth',2);
subplot(2,3,2);
plot(wls,cnAborbance(2,:)','k','LineWidth',2);
subplot(2,3,3);
plot(wls,cnAborbance(3,:)','k','LineWidth',2);

subplot(2,3,4); hold on
plot(wls,log10(srAborbance(1,:)'),'r','LineWidth',4);
title('Carroll-Neitz/Stockman-Rider');
subplot(2,3,5); hold on
plot(wls,log10(srAborbance(2,:)'),'g','LineWidth',4);
title('Carroll-Neitz/Stockman-Rider');
subplot(2,3,6); hold on
plot(wls,log10(srAborbance(3,:)'),'b','LineWidth',4);
title('Carroll-Neitz/Stockman-Rider');
subplot(2,3,4); hold on
plot(wls,log10(cnAborbance(1,:)'),'k','LineWidth',2);
subplot(2,3,5);
plot(wls,log10(cnAborbance(2,:)'),'k','LineWidth',2);
subplot(2,3,6);
plot(wls,log10(cnAborbance(3,:)'),'k','LineWidth',2);

% Baylor
coneParams = DefaultConeParams('cie_baylor');
[T,T_energy,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
    ComputeObserverFundamentals(coneParams,S);
bAborbance = adjIndDiffParams.absorbance;

bFigure = figure; clf; 
subplot(2,3,1); hold on
plot(wls,srAborbance(1,:)','r','LineWidth',4);
title('Baylor/Stockman-Rider');
subplot(2,3,2); hold on
plot(wls,srAborbance(2,:)','g','LineWidth',4);
title('Baylor/Stockman-Rider');
subplot(2,3,3); hold on
plot(wls,srAborbance(3,:)','b','LineWidth',4);
title('Baylor/Stockman-Rider');
subplot(2,3,1); hold on
plot(wls,bAborbance(1,:)','k','LineWidth',2);
subplot(2,3,2);
plot(wls,bAborbance(2,:)','k','LineWidth',2);
subplot(2,3,3);
plot(wls,bAborbance(3,:)','k','LineWidth',2);

subplot(2,3,4); hold on
plot(wls,log10(srAborbance(1,:)'),'r','LineWidth',4);
title('Baylor/Stockman-Rider');
subplot(2,3,5); hold on
plot(wls,log10(srAborbance(2,:)'),'g','LineWidth',4);
title('Baylor/Stockman-Rider');
subplot(2,3,6); hold on
plot(wls,log10(srAborbance(3,:)'),'b','LineWidth',4);
title('Baylor/Stockman-Rider');
subplot(2,3,4); hold on
plot(wls,log10(bAborbance(1,:)'),'k','LineWidth',2);
subplot(2,3,5);
plot(wls,log10(bAborbance(2,:)'),'k','LineWidth',2);
subplot(2,3,6);
plot(wls,log10(bAborbance(3,:)'),'k','LineWidth',2);
%}

