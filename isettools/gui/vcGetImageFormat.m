function iFormat = vcGetImageFormat(data, wave)
% Determine the ISET image format, either RGB or XW, from data.
%
% Syntax:
%   iFormat = vcGetImageFormat(data, wave)
%
% Description:
%    All 3D arrays are RGB format, and the length of wave must equal the
%    3rd dimension.
%
%    The algorithm is this: If DATA is a 3D matrix and the length of wave
%    matches the third dimension, then iFormat is RGB.
%    If DATA is a 2D matrix and length(wave) is 1, then iFormat is RGB.
%
%    If DATA is an NxM matrix and the wave has dimension m, then
%    iFormat'XW' is returned.
%    If DATA is Nx1 and wave has n-entries, then XW is returned.
%
%    The code below contains examples of function usage. To access, type
%    'edit vcGetImageFormat.m' into the Command Window.
%
% Inputs:
%    data    - Matrix. XW or RGB formatted spectral power distribution.
%    wave    - Vector. Wavelength vector.
%
% Outputs:
%    iFormat - String. String containing image format. 'RGB' or 'XW'.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    changeColorSpace (another way to handle some image format issues.)
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    05/09/18  jnm  Formatting

% Examples:
%{
    data = rand(10, 10, 31);
    wave = 400:10:700;
    vcGetImageFormat(data, wave)
    data = RGB2XWFormat(data);
    vcGetImageFormat(data, wave);
%}

if ndims(data) == 3 && length(wave) == size(data, 3)
    % Checks out as classic
    iFormat = 'RGB';
elseif ismatrix(data) && length(wave) == 1
    % A matrix (image) with one wavelength
    iFormat = 'RGB';
elseif ismatrix(data) && (length(wave) == size(data, 2))
    % A matrix as typical XW with columns -- wave
    iFormat = 'XW';
elseif length(data) == length(wave)
    % The data is a vector with the same number of entries as wave. We
    % assume this is XW for one point, and we make sure the data are
    % formatted so things will
    iFormat = 'XW';
else
    iFormat = [];
    s = size(data);
    fprintf('Unrecognized image format. data dimension [%.0f]\n', s)
    fprintf('nWavelength %.0f', length(wave));
end

end
