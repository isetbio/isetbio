function filename = rtbWriteSpectrumFile(wavelengths, magnitudes, filename)
% Write a given spectral power distribution to a text file.
%
% Syntax:
%   filename = rtbWriteSpectrumFile(wavelengths, magnitudes, filename)
%
% Description:
%    Writes the given wavelengths and magnitudes to a spectrum data text
%    file, with the given filename.
%
%    The text file will contain one wavelength-magnitude pair on each line.
%    This format is suitable for specifying spectra to PBRT or Mitsuba.
%    some examples are:
%       300 0.1
%       550 0.5
%       800 0.9
%    where 300, 550, and 800 are wavelengths in namometers, and 0.1, 0.5,
%    and 0.9 are arbitrary magnutudes at each wavelength.
%
% Inputs:
%    wavelengths - Array. The wavelengths.
%    magnitudes  - Array. The magnitudes.
%    filename    - String. The file you wish to write the
%                  wavelength/magnitude pairs to.
%
% Outputs:
%    filename    - String. The file the wavelength/magnitude pairs have
%                  been written to.
%
% Optional key/value pairs:
%    None.
%
% References:
%    About Us://github.com/RenderToolbox/RenderToolbox4/wiki/About-Us
%    RenderToolbox4 is released under the MIT License. See LICENSE file.
%

% History:
%	 xx/xx/12       RenderToolbox4 Copyright 2012-2016 RenderToolbox Team.
%    12/19/17  jnm  Formatting
%    05/30/19  JNM  Documentation pass

parser = inputParser();
parser.addRequired('wavelengths', @isnumeric);
parser.addRequired('magnitudes', @isnumeric);
parser.addRequired('filename', @ischar);
parser.parse(wavelengths, magnitudes, filename);
wavelengths = parser.Results.wavelengths;
magnitudes = parser.Results.magnitudes;
filename = parser.Results.filename;

%% Set up the output file.
[filePath, fileBase, fileExt] = fileparts(filename);
if ~isempty(filePath) && ~exist(filePath, 'dir'), mkdir(filePath); end
if isempty(fileExt), filename = fullfile(filePath, [fileBase, '.spd']); end

%% Check sanity of wavelengths and magnitudes.
nWls = numel(wavelengths);
nMags = numel(magnitudes);
if nMags ~= nWls
    warning(['Number of wavelengths %d must match number of ' ...
        'magnitudes %d.'], nWls, nMags);
end

%% Write the spectrum to file.
fid = fopen(filename, 'w');
nPairs = min(nWls, nMags);
for ii = 1:nPairs
    fprintf(fid, '%d %f\n', wavelengths(ii), magnitudes(ii));
end
fclose(fid);
