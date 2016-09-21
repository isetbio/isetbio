function str = description(obj,varargin)
% Object description for the cone mosaic HEX case
% 
% Adds hex specific information to the base class information

p = inputParser;
p.addRequired('obj', @(x) isa(x, 'coneMosaicHex'));
p.addParameter('skipPigment', false, @islogical);
p.addParameter('skipMacular', false, @islogical);
p.parse(obj, varargin{:});

% String for the base class
str = description@coneMosaic(obj);

% Now add Hex class string
str = addText(str,sprintf('  **Hex info**  \n%s %2.0f\n', 'Resampling factor:', obj.resamplingFactor));
str = addText(str,sprintf('%s %2.0f cols x %2.0f rows\n', 'Rectangular grid:', size(obj.patternOriginatingRectGrid,2), size(obj.patternOriginatingRectGrid,1)));
str = addText(str,sprintf('%s %2.0f cols x %2.0f rows\n', 'Resampled grid:', obj.cols, obj.rows));

end
