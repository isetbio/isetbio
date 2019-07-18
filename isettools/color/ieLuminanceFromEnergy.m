function lum = ieLuminanceFromEnergy(energy, wave, varargin)
% Calculate luminance and related quantities from spectral energy
%
% Syntax:
%   lum = ieLuminanceFromEnergy(energy, wave, [varargin])
%
% Description:
%    Calculate luminance (cd/m2) and related quantities (lux, lumens, cd)
%    from spectral power distribution in energy units.
%
%    The CIE formula for luminance converts a spectral radiance
%    distribution (W/[m2-sr-nm]) into luminance (candelas per meter
%    squared, cd/m2). This routine accepts RGB or XW (space-wavelength)
%    formatted inputs. In XW format, the spectral distributions are in the
%    rows of the ENERGY matrix.
%
%    The formula for luminance and illuminance are the same, differing only
%    in the units of the input. Hence, this routine calculates illuminance
%    (lux) from a spectral irradiance distribution (W/m2-nm). It also
%    calculates luminous intensity (cd) from spectral radiant intensity
%    (W/sr-nm); finally, it calculates luminous flux (lumens, lm) from
%    spectral power (W/nm). The pairings are:
%        Luminance:          cd/m2  from W/[sr-m2-nm]
%        Illuminance:        lux  from  W/[m2-nm]
%        Luminous flux:      lumens from W/nm
%        Luminous intensity: cd from W/[sr-nm]
%
%    To calculate luminance (or illuminance) from a spectral power
%    distribution in photons, use ieLuminanceFromPhotons()
%
%    This function contains examples of usage inline. To access these, type
%    'edit ieLuminanceFromEnergy.m' into the Command Window.
%
% Inputs:
%    energy  - Matrix. The spectral radiance distribution in
%              W / (m^2 * sr * nm), or other units as in description above.
%    wave    - Vector. The wavelengths.
%
% Outputs:
%    lum     - Vector. The luminance, in candelas per meter squared, or
%              other units as in description above.
%
% Optional key/value pairs:
%    'quiet' - Boolean. A boolean indicating whether or not to suppress
%              diagnositic printout to command window. Default is true.
%
% References:
%   http://www.optics.arizona.edu/Palmer/rpfaq/rpfaq.htm
%
% See Also:
%   ieLuminanceFromPhotons
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/30/17  jnm  Comments & formatting
%    11/16/17  jnm  Formatting
%    07/11/19  JNM  Formatting update

% Examples:
%{
   wave = 400:10:700;
   energy = [blackbody(wave, 5000, 'energy'), ...
        blackbody(wave, 6000, 'energy')];
   energy = energy';
   lum = ieLuminanceFromEnergy(energy, wave)
%}

%% Parse input
p = inputParser;
p.addRequired('energy', @isnumeric);
p.addRequired('wave', @isnumeric);
p.addParameter('quiet', true, @islogical);
p.parse(energy, wave, varargin{:});
quiet = p.Results.quiet;

%% Convert data shape
% xwData = ieConvert2XW(energy, wave);
switch vcGetImageFormat(energy, wave)
    case 'RGB'
        xwData = RGB2XWFormat(energy);
    otherwise
        % XW format
        xwData = energy;
end

fName = fullfile(isetbioDataPath, 'human', 'luminosity.mat');
V = ieReadSpectra(fName, wave);

% 683 is the standard factor for conversion when the energy are in Watts.
% The wavelength difference accounts for the wavelength sampling.
if numel(wave) > 1
    dWave = wave(2) - wave(1);
else
    dWave = 10;
    if ~quiet
        disp('ieLuminanceFromEnergy monochrome: 10 nm band assumed');
    end
end
lum = 683 * (xwData * V) * dWave;

end