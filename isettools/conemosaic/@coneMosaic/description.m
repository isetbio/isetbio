function str = description(obj, varargin)
% Create string with key mosaic parameters
%
% We normally include these because they are part of the
% mosaic object. But, we can skip if we like.
%
% Optional parameters (Name-Value pairs)
%   'skipPigment'  - whether to include pigment description
%   'skipMacular'  - whether to include macular description
%
% Output:
%   str - description string

%% parse input
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaic'));
p.addParameter('skipPigment', false, @islogical);
p.addParameter('skipMacular', false, @islogical);
p.parse(obj, varargin{:});

%% Build string using NP's formatting
str = sprintf('%s %2.1f (w) x %2.1f (h)\n', 'Mosaic size (mm):', obj.width*1e3, obj.height*1e3);
txt = sprintf('%s %2.2f (w) x %2.2f (h) \n', 'FOV (deg):', obj.fov(1), obj.fov(2));
str = addText(str,txt);
txt = sprintf('%s %2.2f (w) x %2.2f (h)\n', 'Aperture (um):', obj.pigment.width*1e6, obj.pigment.height*1e6);
str = addText(str,txt);
txt = sprintf('%s %d\n', 'Active cones:' , numel(find(obj.pattern > 1)));
str = addText(str,txt);
txt = sprintf('%s %g\n', 'Density (cones/mm^2):', numel(find(obj.pattern > 1))/(obj.width*obj.height*1e6));
str = addText(str,txt);
txt = sprintf('%s %g\n', 'N Time samples:', length(obj.absorptionsTimeAxis));
str = addText(str,txt);

% cone pigment properties
if ~p.Results.skipPigment
    str = addText(str, obj.pigment.description);
end

% macular pigment properties
if ~p.Results.skipMacular
    str = addText(str, obj.macular.description);
end

end