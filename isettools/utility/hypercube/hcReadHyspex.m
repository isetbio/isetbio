function [img, info] = hcReadHyspex(filename, lines, samples, bands)
% Reads an ENVI image.
%
% Syntax:
%   [img, info] = hcReadHyspex(filename, lines, samples, bands)
%
% Description:
%    Reads the ENVI image in filename and returns the image img as well as
%    the header information in the struct info. Matlab's Multibandread is
%    used to read the data.
%
%    [img, info] = hcReadHyspex(filename, lines, samples, bands) reads only
%    the lines, samples and bands specified in the arguments.
%
%    Renamed for use in ISET-4.0. Taken from read_ENVI_img
%
% Inputs:
%    filename - The string containing the file's name.
%    lines    - (Optional) The total number of rows. Default is [], which
%               is then overwritten by what is in the provided file.
%    samples  - (Optional) The total number of elements in each line.
%               Default is [], which is then overwritten by what is in the
%               provided file.
%    bands    - (Optional) The dimensions of the data when read in its
%               entirety. Default is [], which is then overwritten by what
%               is in the provided file.
%
% Outputs:
%    img      - The hypercube image
%    info     - Associated image header information
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: The uint16 values from the hyspex are not in meaningful energy
%      units (though they are in energy, not photons). Let's see if we can
%      get the scale factor from them that tells us how to scale the uint16
%      values to real units. If we can't, let's do something plausible.
%
% See Also:
%    multiBandRead, hcReadHyspexImginfo.
%

% History:
%    xx/xx/xx  th   Created. Author: trym.haavardsholm@ffi.no
%    12/06/17  jnm  Formatting
%    01/26/18  jnm  Formatting update to match Wiki.

if nargin == 1
    lines = [];
    samples = [];
    bands = [];
elseif nargin == 2
    samples = [];
    bands = [];
elseif nargin == 3
    bands = [];
elseif nargin ~= 4
    error('Wrong number of arguments!');
end

info = hcReadHyspexImginfo(filename);

if isempty(lines) && isempty(samples) && isempty(bands)
    img = multibandread(filename, ...
                        [info.lines, info.samples, info.bands], ...
                        info.data_type, ...
                        info.header_offset, ...
                        info.interleave, ...
                        info.byte_order);
                    
elseif ~isempty(lines) && isempty(samples) && isempty(bands)
    img = multibandread(filename, ...
                        [info.lines, info.samples, info.bands], ...
                        info.data_type, ...
                        info.header_offset, ...
                        info.interleave, ...
                        info.byte_order, ...
                        {'Row', 'Direct', lines});
                    
else
    if isempty(lines), lines = 1:info.lines; end
    if isempty(samples), samples = 1:info.samples; end
    if isempty(bands)
        bands = 1:info.bands;
    elseif strcmpi(bands, 'default')
        bands = info.default_bands;
    end

    img = multibandread(filename, ...
                        [info.lines, info.samples, info.bands], ...
                        info.data_type, ...
                        info.header_offset, ...
                        info.interleave, ...
                        info.byte_order, ...
                        {'Row', 'Direct', lines}, ...
                        {'Column', 'Direct', samples}, ...
                        {'Band', 'Direct', bands});
end