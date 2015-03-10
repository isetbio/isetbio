function coneContrast = humanConeContrast(signalSPD,backgroundSPD,wave,units,mpDensity)
%Calculate  cone contrast of signalSPD superimposed on backgroundSPD
%
%   maxConeContrast = humanConeContrast(signalSPD,backgroundSPD, 
%                                   wave,[units='energy'],[mpDensity=[]])
%
%   The signal and background are in radiance units of energy
%   (watts/sr/m2/nm) or photons (q/sr/m2/nm).  The sample wavelengths are
%   specified in wave (nm). The Stockman cones with a specified macular
%   pigment density mpDensity are used for the calculation.  
%
% Example:
%   coneContrast = humanConeContrast(signalSPD,bgSPD, 400:1:700,'photons')
%   coneContrast = humanConeContrast(signalSPD,bgSPD,380:4:730,'energy')
%   coneContrast = humanConeContrast(signalSPD,bgSPD,380:4:730,'energy',0)
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('signalSPD'),     error('Signal spd required'); end
if notDefined('backgroundSPD'), error('Background spd required'); end
if notDefined('wave'),       error('Wavelength values required'); end
if notDefined('units'),     units = 'energy'; end
if notDefined('mpDensity'), mpDensity = []; end

if length(wave) ~= length(signalSPD), error('Wavelength incorrect.'); end

% Calculation must take place in energy units
if strcmp(units,'photons') || strcmp(units,'quanta')
    signalSPD     = Quanta2Energy(wave,signalSPD(:)');
    backgroundSPD = Quanta2Energy(wave,backgroundSPD');
end

% The max nominal cone contrast is computed here.  We correct for the
% fact that this is less than 100% later.
cones = humanCones('stockmanAbs',wave,mpDensity);
backCones = cones'*backgroundSPD(:);  % The background is always positive
sigCones = cones'*(signalSPD(:)); % Signal can be an increment or decrement
coneContrast = sigCones ./ backCones;

end