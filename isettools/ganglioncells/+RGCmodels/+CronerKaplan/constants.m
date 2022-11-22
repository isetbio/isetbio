classdef constants
% Constants defined in Croner&Kaplan (1995): ''Receptive fields of P and M ganglion cells across the primate retina'

    properties (Constant)
        % Characteristic surround radius as a function of eccentricity (Figure 4 caption)
        % Note: this model was fitted to data of both the M and P RFs
        surroundCharacteristicRadiusFromFitToPandMcells = @(eccDegs) max([0.1*ones(numel(eccDegs),1), 0.203 * eccDegs(:).^0.472],[],2);
        
        % Center peak sensitivity from center characteristic radius in degs
        % (Figure 5b caption)
        centerPeakSensitivityFromCharacteristicRadiusDegsForPcells = @(characteristicRadiusDegs) (0.391 * characteristicRadiusDegs.^(-1.1850));
        
        % Surround peak sensitivity from surround characteristic radius in degs
        % (Figure 5c caption)
        surroundPeakSensitivityFromCharacteristicRadiusDegsForPcells = @(characteristicRadiusDegs) (0.128 * characteristicRadiusDegs.^(-2.147));
        
        % Surround / center Rc ratio
        surroundToCenterRcRatio = 6.7;

        % From Table 1, Pcells (0-40 degs)
        surroundToCenterPeakSensitivityRatio = 1.1/106.3;

        % Surround radius from center radius
        surroundRadiusFromCenterRadiusDegsForPcells = @(centerRadiusDegs) (centerRadiusDegs * RGCmodels.CronerKaplan.constants.surroundToCenterRcRatio);
        
        % Surround-to-Center integrated sensitivity ratio (Figure 11 caption)
        surroundToCenterIntegratedSensitivityRatioFromEccDegsForPcells = @(eccDegs)  (0.466 + 0.007 * abs(eccDegs));
    end
end