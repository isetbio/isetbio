function str = description(obj, varargin)
% Object description for the cone mosaic HEX case
%
% Syntax:
%   str = description(obj, varargin)
%
% Description:
%    Create the object description for the cone mosaic hex object.
%
%    Adds hex specific information to the base class information
%
%    See arguments in the base class for the pigment properties such as:
%    skipMacularPigment, and skipPhotopigment
%
% Inputs:
%    obj      - The cone mosaic hex object
%    varargin - (Optional) Any additional information to add. Default empty
%
% Outputs:
%    str      - The string containing the cone mosaic hex's description.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/16  BW   ISETBIO Team, 2016
%    02/16/18  jnm  Formatting

%% Parse
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaicHex'));

%% String for the base class.
% Argument parsing happens there.
str = description@coneMosaic(obj, varargin{:});

%% Now add Hex class string
str = addText(str, sprintf('  **Hex info**  \n%s %2.0f\n', ...
    'Resampling factor:', obj.resamplingFactor));
str = addText(str, sprintf('%s %2.0f cols x %2.0f rows\n', ...
    'Rectangular grid:', size(obj.patternOriginatingRectGrid, 2), ...
    size(obj.patternOriginatingRectGrid, 1)));
str = addText(str, sprintf('%s %2.0f cols x %2.0f rows\n', ...
    'Resampled grid:', obj.cols, obj.rows));

end
