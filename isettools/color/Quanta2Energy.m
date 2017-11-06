function energy = Quanta2Energy(wavelength, photons)
% Convert quanta (photons) to energy (watts)
%
% Syntax:
%   energy = Quanta2Energy(wavelength, photons)
%
% Description:
%    Convert photons represented at the sampled wavelength positions to
%    energy units (watts or joules).  You get watts if photons specifies
%    photon rate (photons/sec) and joules if it just specifies the number
%    of photons.
%
% Inputs:
%    wavelength - Column vector describing the wavelength samples [nm]
%    photons    - Matrix in either RGB or XW (space-wavelength) format.
%                 In the XW format each spatial position is in a row and
%                 the wavelength varies across columsn. The output, ENERGY, 
%                 [watts or Joules] is returned in  same format as input
%                 (RGB or XW). CAUTION: The input form differs from the 
%                 Energy2Quanta() call, which uses the energy spectra in
%                 the columns.
%
% Outputs:
%    energy     - The energy, in watts or joules, is returned in the same
%                 format that the photon matrix was in (RGB or XW)
%
% Notes:
%    * [NOTE: XXX - CAUTION: The input form differs from the
%      Energy2Quanta() call, which has the energy spectra in the columns.]
%    * [NOTE: XXX - We should regularize the calls to Energy2Quanta() and
%      this routine, probably by making the other routine take RGB or XW
%      format as well. Old legacy issues, sigh.]
%    * [NOTE: DHB - When this is passed a single column vector, it doesn't
%      complain (because that's an illegal format) but rather returns a
%      single row vector.  Not sure whether we should check and throw an
%      error for this case.]
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/27/17  jnm  Comments & formatting
   

% Examples:
%{
   wave = 400:10:700;  
   p = blackbody(wave, 3000:1000:8000, 'photons');
   e = Quanta2Energy(wave, p'); e = diag(1./e(:, 11))*e;
   figure; plot(wave, e')

   p1 = blackbody(wave, 5000, 'photons');
   e = Quanta2Energy(wave, p1');           % e is a row vector in XW format
   p2 = Energy2Quanta(wave, transpose(e)); % Notice the TRANSPOSE
   figure; plot(wave, p1, 'ro', wave, p2, 'k-')
%}

if isempty(photons)
    energy = [];
    return;
end

% Make wavelenth as row vector
wavelength = wavelength(:)'; 

% Fundamental constants
h = vcConstants('h');		% Planck's constant [J sec]
c = vcConstants('c');		% speed of light [m/sec]

% Main routine handles RGB or XW formats
iFormat = vcGetImageFormat(photons, wavelength);

switch iFormat
    case 'RGB'
        [n, m, w] = size(photons);
        if w ~= length(wavelength)
            error('Quanta2Energy: photons third dimension must be nWave');
        end
        photons = RGB2XWFormat(photons);
        
        % energy = (h * c / (1e-9)) * (photons ...
        %     ./ repmat(wavelength, n * m, 1));
        energy = (h * c / 1e-9) * bsxfun(@rdivide, photons, wavelength);
        energy = XW2RGBFormat(energy, n, m);

    case 'XW'
        % If photons is a vector, it must be a row
        if isvector(photons)
            photons = photons(:)';
        end
        if size(photons, 2) ~= length(wavelength)
            error('Quanta2Energy: quanta must have length of nWave');
        end
        energy = (h * c / 1e-9) * bsxfun(@rdivide, photons, wavelength);
        
    otherwise
        error('Unknown image format');
end
end
