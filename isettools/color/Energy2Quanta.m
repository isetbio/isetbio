function photons = Energy2Quanta(wavelength, energy)
% Convert energy (watts or joules) to photon units
%
% Syntax:
%   photons = Energy2Quanta(wavelength, energy)
%
% Description
%    Energy (e.g., watts) is returned in units of photons, that is from
%    watts/sr/m2/nm to photons/sec/sr/m2/nm. The parameter wavelength is
%    the vector of sample wavelengths in units of nanometers (nm). If you
%    pass watts, you get photons/sec, if you pass joules you get photons
%
%    The return (photons) is in the same format.
%
%    The columns of the matrix energy are different spatial samples. The
%    entries down a column represent the energy as a function of
%    wavelength. Thus, the entries across a row represent the energy at a
%    single wavelength for each of the samples. This format is unfortunate,
%    because it is the transpose of the standard XW format used throughout
%    the rest of ISET. It is also the transpose of the format used in the
%    function Quanta2Energy.
%
%    This function contains examples of usage inline. To access these, type
%    'edit Energy2Quanta.m' into the Command Window.
%
% Inputs:
%    wavelength - Vector. A vector of sample wavelengths in nm.
%    energy     - Matrix. A matrix containing the spatial samples. Each
%                 column represents one spatial sample, with each row
%                 corresponding to one wavelength (transpose of XW format).
%
% Outputs:
%    photons    - Matrix. A matrix in the same format as input (transpose
%                 of XW format).
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [NOTE - BAW: In the fullness of time, this function will be
%      converted to XW format; but it will require some effort because they
%      are inserted in so many important places. I am going to need a long
%      holiday to make this change. Currently, this function is not
%      transparently compatible with its counterpart, Quanta2Energy. If
%      this ever happens, fix comments above and in the example.]
%
% See Also:
%   Quanta2Energy
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/30/17  jnm  Comments & formatting
%    11/16/17  jnm  Formatting
%    07/11/19  JNM  Formatting update

% Examples:
%{
    wave = 400:10:700;
    in = [blackbody(wave, 5000, 'energy'), ...
        blackbody(wave, 6000, 'energy')];
    p = Energy2Quanta(wave, in);
    figure;
    plot(wave, p);

    % Notice that we have to transpose p when we call the inverse
    % Quanta2Energy. Tragic.

    % Now, from Quanta2Energy, out in XW format
    out = Quanta2Energy(wave, p');
    figure;
    plot(wave, in, 'ro', wave, out', 'k-');
%}

if ~isvector(wavelength)
    error('Energy2Quanta: Wavelength must be a vector');
else
    wavelength = wavelength(:);  % Force to column
end

% Fundamental constants.
h = vcConstants('h');  % Planck's constant [J sec]
c = vcConstants('c');  % speed of light [m/sec]

if ndims(energy) == 3
    [n, m, w] = size(energy);
    if w ~= length(wavelength)
        error('%s:  size(energy, 3) must equal numWave', mfilename);
    end
    energy = reshape(energy, n * m, w)';
    photons = (energy / (h * c)) .* (1e-9 * wavelength(:, ones(1, n * m)));
    photons = reshape(photons', n, m, w);
else
    [n, m] = size(energy);
    if n ~= length(wavelength)
        error('%s: energy must have length of numWave', mfilename);
    end
    photons = (energy / (h * c)) .* (1e-9 * wavelength(:, ones(1, m)));
end

end