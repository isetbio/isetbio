function isomerizations = IsomerizationsFromLuminanceGeisler(luminanceCdM2,durationInSeconds,pupilDiameterMm,varargin)
% Estimate cone isomerizations from luminance level
%
% Syntax:
%    isomerizations = ...
%     IsomerizationsFromLuminanceGeisler(luminanceCdM2, ...
%                 durationInSeconds,pupilDiameterMm)
%
% Description:
%    Geisler 1(984) gives a formula for isomerizations from luminance
%    (Equation 2, p. 776), obtained from Wyszecki and Stiles.
% 
%    These estimates are reasonable for L and M cones, but not for S
%    cones.
%
%    This function computes using the 1984 parameters, ignoring convolution
%    by the optical point spread function but including lens transmittance.
%
% References:
%    Geisler, W.S. 1984 Physical limits of acuity and hyperacuity. Journal
%    of the Optical Society of America A 1, 775-782.
%
% Inputs:
%    luminanceCdM2      - Luminance in cd/m2
%    durationInSeconds  - Duration over which to compute isomerizations
%    pupilDiameterMm    - Pupil diameter in mm
%
% Outputs:
%    isomerizations     - Estimate of number of isomerizations.
%
% Optional key/value pairs:
%    'coneApertureDiameterMinutes' - Collection aperture of a cone in minutes
%                                    of arc. Default 0.6.
%    'coneQuantalEfficiency555'    - Cone quantal efficiency at 555 nm.
%                                    Default 0.5;
%    'occularTransmittance'        - Fraction of light transmitted through
%                                    preretinal absorption. Default 0.68.
%
% See also: v_IrradianceIsomerizations

% History:
%   05/12/18  dhb  Brought over from IBIOColorDetect.

% Examples:
%{
   % 5 cd/m2 -> 1000 isomerizations/sec seems ballpark right.
   lum = 5;  % cd/m2
   dur = 1;  % sec
   pupilDiameter = 3; % mm
   isomerizations = IsomerizationsFromLuminanceGeisler(lum,dur,pupilDiameter)
%}

%% Parse input
p = inputParser;
p.addParameter('coneApertureDiameterMinutes',0.6,@isnumeric);
p.addParameter('coneQuantalEfficiency555',0.5,@isnumeric);
p.addParameter('occularTransmittance',0.68,@isnumeric);
p.parse(varargin{:});

%% Geisler's choice of occular transimttance parameter
occularTransmittance = p.Results.occularTransmittance;

% Geisler's choice of quantum efficiency at 555 nm parameter
quantumEfficiency555 = p.Results.coneQuantalEfficiency555;

% Cone acceptance apperture
coneApertureDiamterMinutes = p.Results.coneApertureDiameterMinutes;
coneAreaMinutes2 = pi*((coneApertureDiamterMinutes/2)^2);

% Pupil area
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);

% This is a magic unit constant in the formula
magicUnitConstant = 347.8;

% Apply the formula
isomerizations = coneAreaMinutes2*durationInSeconds* ...
    pupilAreaMm2*occularTransmittance* ...
    quantumEfficiency555*magicUnitConstant*luminanceCdM2;

end