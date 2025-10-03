function [Kwide, Knarrow, Rnarrow] = doubleExponentSurroundRFparams(...
    Kc, Rwide, KsToKcPeakRatio, narrowToWideVolumeRatio, RnarrowToRwideRatio)
% Given Kc (center gain), KsToKcRatio (surround/center peak sensitivity
% ratio) and Rwide (wide RF surround radius), narrowToWideVolumeRatio and RnarrowToRwideRatio,
% compute the peak sensitivity gains for the wide and the narrow surround components

    % Params for the sum of 2 exponentials (wide field + narrow field, Packer& Dacey 2002)
    % RF = Kwide * exp(-2.3*R/Rwide) + Knarrow * exp(-2.3*R/Rnarrow) 
    % Rwide: radius at which sensitivity drops to 10%, which is defined as half the RFdiameter
    % Volume: R^2 * K

    Kwide = KsToKcPeakRatio*Kc / (1 + narrowToWideVolumeRatio / (RnarrowToRwideRatio^2));
    Knarrow  = Kwide * narrowToWideVolumeRatio / (RnarrowToRwideRatio^2);
    Rnarrow  = RnarrowToRwideRatio * Rwide;
end