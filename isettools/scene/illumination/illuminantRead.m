function spectralRadiance = illuminantRead(illP, lightName)
% Return spectral radiance of a standard illuminants in energy units
%
% Syntax:
%	spectralRadiance = illuminantRead([illP], [lightName])
% 
% Description:
%    The illuminant parameters (illP) are stored in the structure illP. See
%    the example below for how to initialize the values of this structure.
%    The illP structure is idiosyncratic and used only here. 
% 
%    If you don't wish to establish illP, but only to get a default SPD for
%    some named illuminant, you can use the format:
%        illuminantRead([], 'd65')
%
%    In this case the spectral radiance is returned at 400:10:700 nm
%    samples and the mean luminance is 100 cd/m2.
%
%    There are examples contained in the code. To access the examples, type
%    'edit illuminantRead.m' into the Command Window.
%
%    The standard illuminant names are:
%       {'tungsten'}
%       {'illuminantc'}
%       {'d50'}
%       {'fluorescent'}
%       {'d65', 'D65'}
%       {'equalenergy'}
%       {'blackbody'}   -- You must specify a color temperature in
%                          illP.temperature
%       {'555nm'}
%
% Inputs:
%    illP             - (Optional) The illuminant parameters structure.
%                       Contains values such as: name, spectrum.wave,
%                       luminance, temperature. Default is providing the
%                       following information:
%                           wave = 400:10:700
%                           luminance = 100
%    lightName        - (Optional) The name of the illuminant. Options are:
%                       555nm, blackbody*, d50, d65, equal energy, equal
%                       photons, fluorescent, illuminant c, and tungsten.
%                       Illuminants with (*) require their color
%                       temperature be specified with illP.temperature.
%                       Default is 'd65'.
%
% Outputs:
%    spectralRadiance - The requested spectral radiance
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    illuminantCreate, illuminantGet/Set
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%	 01/30/18  jnm  Formatting

% Examples:
%{
    illuminantRead([], 'd65')
    illSPD = illuminantRead([], 'tungsten');
    plot(400:10:700, illSPD)
%}
%{
    illP.name = 'd65';
    illP.spectrum.wave = 400:10:700;
    illP.luminance = 100;
    plot(illuminantRead(illP));
%}
%{
    illP.name = 'blackbody';
    illP.temperature = 3000;
    illP.spectrum.wave = 400:10:700;
    illP.luminance = 100;
    sr = illuminantRead(illP);
    plot(sr)
%}

if notDefined('illP')
    if notDefined('lightName')
        name = 'd65';
        warndlg('No illuminant name. Assuming D65');
    else
        name =  lightName;
    end
    luminance = 100;
    wave = 400:10:700;
else
    name = illP.name;
    luminance = illP.luminance; 
    wave = illP.spectrum.wave;
end

name = ieParamFormat(name);

switch name
    case {'tungsten'}
        SPD = ieReadSpectra('data/lights/Tungsten', wave);
    case {'illuminantc'}
        SPD = ieReadSpectra('data/lights/illuminantC', wave);
    case {'d50'}
        SPD = ieReadSpectra('data/lights/D50', wave);
    case {'fluorescent'}
        SPD = ieReadSpectra('data/lights/Fluorescent', wave);
    case {'d65'}
        SPD = ieReadSpectra('data/lights/D65', wave);
    case {'white', 'uniform', 'equalenergy'}
        SPD = ones(length(wave), 1);
    case {'equalphotons'}
        SPD = Quanta2Energy(wave, ones(1, length(wave)))';
    case 'blackbody'
        if ~checkfields(illP, 'temperature')
            temperature = 6500;
        else
            temperature = illP.temperature;
        end
        SPD = blackbody(wave, temperature);
    case {'555nm', 'monochrome'}
        SPD = zeros(length(wave), 1);
        % Set the wavelength closest to 555 to 1
        [~, idx] = min(abs(wave - 555));
        SPD(idx) = 1;
    otherwise   
        error('Illumination:  Unknown light source');
end

% Compute the current light source luminance; scale it to the desired
% luminance. The formula for luminance is:
% currentL = 683 * binwidth * (photopicLuminosity' * SPD);
currentL = ieLuminanceFromEnergy(SPD', wave);
spectralRadiance = (SPD / currentL) * luminance;

% Just check the values
%  ieLuminanceFromEnergy(spectralRadiance', wave)
%  ieXYZFromEnergy(spectralRadiance', wave)

end