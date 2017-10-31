function photons = Energy2Quanta(wavelength, energy)
% Convert energy (watts) to number of photons.
%
% Syntax:
%   photons = Energy2Quanta(wavelength, energy)
%
% Description
%    Energy (e.g., watts) is returned in units of photons, that is from
%    watts/sr/m2/nm to photons/sr/m2/nm.  The parameter wavelength is the
%    vector of sample wavelengths in units of nanometers (nm).
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
%
% Inputs:
%    wavelength - the vector of sample wavelengths in nm.
%    energy     - A matrix containing the spatial samples. Rows represent
%                 energy at a single wavelength, and columns as a function
%                 of the wavelengths. (Transpose of XW Format)
%
% Outputs:
%    photons    - Matrix in either RGB or XW (space-wavelength) format.
%                 In the XW format each spatial position is in a row and
%                 the wavelength varies across columsn. 
%
% Notes:
%    - XXX: In the fullness of time, this function will be converted to XW
%      format; but it will require some effort because they are inserted in
%      so many important places.  I am going to need a long holiday to make
%      this change.
%
% See Also:
%   Quanta2Energy
%

% History: 
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/30/17  jnm  Comments & formatting

% Examples:
%{
   wave = 400:10:700;  
   in = [blackbody(wave, 5000, 'energy'), blackbody(wave, 6000, 'energy')];
   p = Energy2Quanta(wave, in);  
   figure; plot(wave, p);

   % Notice that in the return, out becomes a row vector, consistent with
   % XW format. Also, notice that we have to transpose p to make this work.
   % Tragic.

   out = Quanta2Energy(wave, p');      % out is a row vector, XW format
   figure; plot(wave, in, 'ro', wave, out', 'k-') 
%}

if ~isvector(wavelength)
        error('Energy2Quanta: Wavelength must be a vector');
else
    wavelength = wavelength(:);      % Force to column
end

% Fundamental constants.  
h = vcConstants('h');		% Planck's constant [J sec]
c = vcConstants('c');		% speed of light [m/sec]

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