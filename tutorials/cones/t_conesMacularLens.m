% Illustrate cone quantal efficiency changes with macular pigment and lens
%
% Description:
%   This tutorial illustrates how to use machinery built into the @cMosaic
%   object to calculate how the cone quantal efficiency changes with
%   eccentricity due to changes in MP density with eccentricity. The effect
%   of the lens on cone quantal efficiency is also illustrated.
%

% History:
%    01/16/23  NPC  ISETBIO Team, Copyright 2023 Wrote it.

function t_MPandLensAdustedQuantalEfficiency

    % Generate a @cMosaic at some desired eccentricity
    targetEccXYdegs = [-1 1];

    sigmaGaussian = 0.204;
    coneApertureModifiers = struct(...
        'smoothLocalVariations', true, ...
        'sigma',  sigmaGaussian, ...
        'shape', 'Gaussian');

    % We will test at (-2.0, 0.0) degs
    theConeMosaic = cMosaic(...
                'eccentricityDegs', targetEccXYdegs, ...
                'sizeDegs', [1 0.5], ...
                'rodIntrusionAdjustedConeAperture', true, ...
                'coneApertureModifiers', coneApertureModifiers);
      

    % Obtain the macular pigment boost factors that must be applied to the
    % retinal image so as to yield an effective cone quantal efficiency
    % that reflects the change in macular pigment with eccentricity
    %
    if (theConeMosaic.eccVaryingMacularPigmentDensity)
        % Each cone has a different boost factor, based on its eccentricity
        macularPigmentBoostFactors = cMosaic.macularPigmentBoostFactors(theConeMosaic.macular, theConeMosaic.coneRFpositionsDegs);
    else
        % Single factor based on the eccentricity of the mosaic
        macularPigmentBoostFactors = cMosaic.macularPigmentBoostFactors(theConeMosaic.macular, theConeMosaic.eccentricityDegs);
    end

    % Plot the distribution of these boost factors (for the peak
    % wavelength) across the cone mosaic
    plotMPBoostMap(theConeMosaic, macularPigmentBoostFactors);

    % Obtain the boost vector at the center of the mosaic
    targetEcc = mean(theConeMosaic.coneRFpositionsDegs,1);
    [~, targetEccIndex] = min(sum((bsxfun(@minus, theConeMosaic.coneRFpositionsDegs, targetEcc)).^2,2));
    targetEcc = theConeMosaic.coneRFpositionsDegs(targetEccIndex,:);
    boostVector = macularPigmentBoostFactors(targetEccIndex,:);
   

    % Quantal efficiencies with the foveal MP density and no lens
    coneQuantalEfficienciesFovealMP = theConeMosaic.qe;

    % Quantal efficiencies with the MP density at the target ecc and no lens
    coneQuantalEfficienciesEccBasedMP = diag(boostVector) * coneQuantalEfficienciesFovealMP;


    % Quantal efficiencies including the Lens
    theLens = Lens('wave',theConeMosaic.wave);
    lensTransmittance = theLens.transmittance;
    coneQuantalEfficienciesFovealMPwithLens = diag(lensTransmittance)*coneQuantalEfficienciesFovealMP;
    coneQuantalEfficienciesFinal = diag(lensTransmittance)*coneQuantalEfficienciesEccBasedMP;

    % Show the various changes induced in the quantal cone efficiencies
    % by the macular pigment and the lens
    hFig = figure(1); clf
    set(hFig, 'Position', [200 300 800 650], 'Color', [1 1 1]);


    subplot(2,2,1)
    plot(theConeMosaic.wave, coneQuantalEfficienciesFovealMP(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(theConeMosaic.wave, coneQuantalEfficienciesFovealMP(:,2), 'g-', 'LineWidth', 1.5);
    plot(theConeMosaic.wave, coneQuantalEfficienciesFovealMP(:,3), 'b-', 'LineWidth', 1.5);
    set(gca, 'YLim', [0 0.5], 'FontSize', 16);
    grid on
    title('foveal MP without lens')

    subplot(2,2,2)
    plot(theConeMosaic.wave, coneQuantalEfficienciesEccBasedMP(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(theConeMosaic.wave, coneQuantalEfficienciesEccBasedMP(:,2), 'g-', 'LineWidth', 1.5);
    plot(theConeMosaic.wave, coneQuantalEfficienciesEccBasedMP(:,3), 'b-', 'LineWidth', 1.5);
    set(gca, 'YLim', [0 0.5], 'FontSize', 16);
    grid on
    title(sprintf('(%2.3f,%2.3f) degs MP without lens', targetEcc(1), targetEcc(2)));

    subplot(2,2,3)
    plot(theConeMosaic.wave, coneQuantalEfficienciesFovealMPwithLens(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(theConeMosaic.wave, coneQuantalEfficienciesFovealMPwithLens(:,2), 'g-', 'LineWidth', 1.5);
    plot(theConeMosaic.wave, coneQuantalEfficienciesFovealMPwithLens(:,3), 'b-', 'LineWidth', 1.5);
    set(gca, 'YLim', [0 0.5], 'FontSize', 16);
    grid on
    title(sprintf('foveal MP with lens'));

    subplot(2,2,4)
    plot(theConeMosaic.wave, coneQuantalEfficienciesFinal(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(theConeMosaic.wave, coneQuantalEfficienciesFinal(:,2), 'g-', 'LineWidth', 1.5);
    plot(theConeMosaic.wave, coneQuantalEfficienciesFinal(:,3), 'b-', 'LineWidth', 1.5);
    set(gca, 'YLim', [0 0.5], 'FontSize', 16);
    grid on
    title(sprintf('(%2.1f,%2.1f) degs MP with lens', targetEcc(1), targetEcc(2)));



end

function plotMPBoostMap(theConeMosaic, macularPigmentBoostFactors)
    hFig = figure(2); clf;
    set(hFig, 'Position', [700 300 800 850], 'Color', [1 1 1]);

    ax = subplot(2,1,1);
    theConeMosaic.visualize(...
        'figureHandle', hFig, ...
        'axesHandle', ax, ...
        'plotTitle', 'the cone mosaic');

    [~,idx] = max(macularPigmentBoostFactors(:));
    [~, wavelengthIndexOfMaxBoost] = ind2sub(size(macularPigmentBoostFactors), idx);

    plotTitle = sprintf('MP-density adjusted boost factor at %d nm', theConeMosaic.wave(wavelengthIndexOfMaxBoost));
    mpBoostMap = macularPigmentBoostFactors(:, wavelengthIndexOfMaxBoost);
    minBoostFactor = min(mpBoostMap(:));
    maxBoostFactor = max(mpBoostMap(:));

    ax = subplot(2,1,2);
    theConeMosaic.visualize('figureHandle', hFig, 'axesHandle', ax, ...
        'activation', mpBoostMap, ...
        'activationRange', [minBoostFactor maxBoostFactor], ...
        'verticalActivationColorBarInside', true, ...
        'backgroundColor', [0.1 0.1 0.1], ...
        'plotTitle', plotTitle);
   
end