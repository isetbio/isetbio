classdef constants
% Constants defined in Croner&Kaplan (1995): ''Receptive fields of P and M ganglion cells across the primate retina'

    properties (Constant)
        % Characteristic surround radius as a function of eccentricity (Figure 4 caption)
        % Note: this model was fitted to data of both the M and P RFs
        surroundCharacteristicRadiusFromFitToPandMcells = @(eccDegs) max([0.1 0.203 * eccDegs.^0.472]);
    end
end

