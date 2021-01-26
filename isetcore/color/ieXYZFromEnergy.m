function XYZ = ieXYZFromEnergy(energy, wave)
% CIE XYZ values from spectral radiance(w/nm/sr/m2) or irradiance(w/nm/m2)
%
% Syntax:
%   XYZ = ieXYZFromEnergy(energy, wave)
%
% Description:
%    Calculate CIE XYZ values of the spectral radiance or irradiance
%    functions in the variable energy. The input format of energy can be
%    either XW (space-wavelength) or RGB. The wavelength samples of energy
%    are stored in the variable wave.
%
%    For XW format input, the SPDs are in the rows, and the returned XYZ
%    values are in the corresponding rows of the returned matrix. That is,
%    the return comes back in XW format if the input is in XW format.
%
%    The units of Y are candelas/meter-squared if energy is radiance and
%    lux if energy is irradiance, and if the units of radiance are
%    watts/[sr-m2-nm] or the units of irradiance are Watts/[m2-nm].
%
%    This function contains examples of usage inline. To access these, type
%    'edit ieXYZFromEnergy.m' into the Command Window.
%
% Inputs:
%    energy - Matrix. Energy in XW (space-wavelength) or RGB formats.
%    wave   - Vector. Wavelength samples of energy.
%
% Outputs:
%    XYZ    - Matrix. CIE XYZ values in format matched to the input. Y is
%             either in candelas per meter squared (radiance), or in lux
%             (irradiance). Units for radiance are w / [sr * m^2 * nm]),
%             and for irradiance w / [m^2 * nm].
%
% Optional key/value pairs:
%    None.
%
% See Also:
%   ieXYZFromPhotons
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    xx/xx/xx  xxx  Return output in RGB format if input came in that way.
%    10/27/17  jnm  Comments & formatting
%    11/01/17  dhb  Remove comments that RGB return is new, since it isn't
%                   anymore.
%    07/10/19  JNM  Formatting update

% Examples:
%{
    wave = 400:10:700;
    tmp = load('CRT-Dell');
    dsp = tmp.d;
    energy = displayGet(dsp, 'spd', wave);
    energy = energy';
    displayXYZ = ieXYZFromEnergy(energy, wave)

    patchSize = 1;
    macbethChart = sceneCreate('macbeth', patchSize);
    p = sceneGet(macbethChart, 'photons');
    wave = sceneGet(macbethChart, 'wave');
    e = Quanta2Energy(wave, p);
    XYZ = ieXYZFromEnergy(e, wave);
%}

%% Force data into XW format.
if ndims(energy) == 3
    if length(wave) ~= size(energy, 3)
        error('Bad format for input variable energy.');
    end
elseif isvector(energy)
    energy = energy(:)';  % Force to row vector
end

iFormat = vcGetImageFormat(energy, wave);
switch iFormat
    case 'RGB'
        % [rows, cols, w] = size(data);
        [xwData, r, c] = RGB2XWFormat(energy);
        % disp('RGB return')
    otherwise
        % XW format
        xwData = energy;
end

% xwData = ieConvert2XW(energy, wave);
if size(xwData, 2) ~= length(wave)
    error('Problem converting input variable energy into XW format.');
end

% The spectra of the energy points are in the rows of xwData. We load the
% XYZ color matching functions into the columns of S.
S = ieReadSpectra('XYZ', wave);
if numel(wave) > 1
    dWave = wave(2) - wave(1);
else
    dWave = 10;
    disp('10 nm band assumed');
end

% The return value has three columns, [X, Y, Z].
XYZ = 683 * dWave * (xwData * S);

% If it was sent in RGB, return it in RGB
switch iFormat
    case 'RGB'
        XYZ = XW2RGBFormat(XYZ, r, c);
    otherwise
        % XW format, do nothing
end

end
