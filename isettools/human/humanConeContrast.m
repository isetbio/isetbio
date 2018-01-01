function coneContrast = humanConeContrast(signalSPD,backgroundSPD,wave,units)
%Calculate  cone contrast of signalSPD superimposed on backgroundSPD
%
% Syntax
%   maxConeContrast = humanConeContrast(signalSPD,backgroundSPD, 
%                                   wave,[units='energy'])
% Description
%   The signal and background radiance spectral power distributions can be
%   in units of energy (watts/sr/m2/nm) or photons (q/sr/m2/nm).  The
%   sample wavelengths are specified in wave (nm). The Stockman cones with
%   default macular pigment density are used for the calculation.
% 
% Inputs
%   signalSPD
%   backgroundSPD
%   wave
%
% Optional key/val pairs
%   None
%
% Return
%   coneContrast - L,M,S
%
% Note:
%   * We should add an option to calculate with different macular pigment
%   densities.  Should be straightforward, but ..
%
% Copyright ImagEval Consultants, LLC, 2005.

% Example:
%{
   dsp   = displayCreate;
   wave  = displayGet(dsp,'wave');
   bgSPD = displayGet(dsp,'spd')*0.5*ones(3,1);
   [~,sigSPD]   = humanConeIsolating(dsp);
   coneContrast = humanConeContrast(sigSPD,bgSPD, wave,'energy')
%}

%% Deal with parameters
if notDefined('signalSPD'),     error('Signal spd required'); end
if notDefined('backgroundSPD'), error('Background spd required'); end
if notDefined('wave'),          error('Wavelength values required'); end
if notDefined('units'),         units = 'energy'; end

if length(wave) ~= length(signalSPD), error('Wavelength incorrect.'); end

%% Calculation takes place in energy units
if strcmp(units,'photons') || strcmp(units,'quanta')
    signalSPD     = Quanta2Energy(wave,signalSPD(:)');
    backgroundSPD = Quanta2Energy(wave,backgroundSPD');
end

%% Compute the cone contrast

% The contrast is the ratio of the background cone absorptions to the
% signal cone absorptions
coneFile = fullfile(isetbioDataPath,'human','stockman');
cones = ieReadSpectra(coneFile,wave);   % plot(wave,spCones)
        
backCones    = cones'* backgroundSPD(:);  % The background is always positive
sigCones     = cones'* signalSPD;         % Signal can be an increment or decrement

coneContrast = diag(1./backCones)*sigCones;

end