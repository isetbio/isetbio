function str = description(obj, varargin)
% Create string with key mosaic parameters
%
% Syntax:
%   str = description(obj, varargin)
%
% Description:
%    Create a string containing key cone mosaic parameters.
%
% Inputs:
%    obj           - cone mosaic object
%
% Outputs:
%    str           - string describing some things about the mosaic.
%
% Optional key/value pairs:
%    'skipPigment' - Omit info about the photopigment in the string.
%                    Default is false.
%    'skipMacular' - Omit info about the maculular pigment in string.
%                    Default is false.
%
% See Also:
%    coneMosaic, photoPigment, macular.
%

%% parse input
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addParameter('skipPigment', false, @islogical);
p.addParameter('skipMacular', false, @islogical);
p.parse(obj, varargin{:});

%% Build string using NP's formatting
str = sprintf('%s %2.1f (w) x %2.1f (h)\n', 'Mosaic size (mm):', ...
    obj.width * 1e3, obj.height * 1e3);
txt = sprintf('%s %2.2f (w) x %2.2f (h) \n', 'FOV (deg):', obj.fov(1), ...
    obj.fov(2));
str = addText(str, txt);
txt = sprintf('%s %2.2f (w) x %2.2f (h)\n', 'Aperture (um):', ...
    obj.pigment.width * 1e6, obj.pigment.height * 1e6);
str = addText(str, txt);
txt = sprintf('%s %d\n', 'Active cones:' , numel(find(obj.pattern > 1)));
str = addText(str, txt);
txt = sprintf('%s %g\n', 'Density (cones/mm^2):', ...
    numel(find(obj.pattern > 1)) / (obj.width * obj.height * 1e6));
str = addText(str, txt);

% Does the number of absorptions always per integration time? I think so.
% But there is the question of absorption times in the os, too.
nSamples = size(obj.absorptions, 3);   % Number of time samples
T = nSamples * obj.integrationTime;
txt = sprintf('%s %g sec (%d samps)\n', 'Duration: ', T, nSamples);
str = addText(str, txt);

% cone pigment properties
if ~p.Results.skipPigment, str = addText(str, obj.pigment.description); end

% macular pigment properties
if ~p.Results.skipMacular, str = addText(str, obj.macular.description); end

end