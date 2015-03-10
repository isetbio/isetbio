function [CFAletters,CFAnumbers,mp] = sensorDetermineCFA(sensor)
% Determine the CFA organization for an image sensor
%
%   [CFAletters,CFAnumbers,mp] = sensorDetermineCFA(sensor)
%
% - Major change in May, 2010 -
%
% The CFAnumbers is an integer matrix with the row/col size as the sensor.
% The integers in the matrix refer to which color filter (as defined by the
% column number occupied by the color filter in the sensor data structure
% (sensorGet(sensor,'filterSpectra')).  These are the same numbers that are
% used by the field 'pattern' to define the color filter array spatial
% layout.
%
% The CFAletters is an array of the same size, but the CFAnumbers are
% replaced by a string that characterizes the filter appearance. If the
% filterName first letter that is within the list defined in
% sensorColorOrder, then we use that letter.  Otherwise, we assign the
% color k (black). This may encourage better notation.  Or, we may get
% annoyed with this and write a routine to determine a reasonable color
% from the filterSpectra curve.
%
% The return mp contains the RGB values associated with the letters that
% are present in the CFAletters.  If there are three letters, say r g and
% b, then the mp values will be 
%
%        1     0     0
%        0     1     0
%        0     0     1
%
% The set of valid letter names for color filters (and their numerical
% order) is set in the function sensorColorOrder.  Type
%
%   sensorColorOrder
%
% to see the list of letters. In May 2010 the list was:
%
% {'r' 'g' 'b' 'c' 'y' 'm' 'w' 'u' 'x' 'z' 'o' 'k'}
%
% The associated RGB values for each letter are also stored in that
% routine. This mapping can be found from
%
%    [letters, mp] = sensorColorOrder;
%
% Examples:
%   sensor = vcGetObject('sensor');
%   [cfa,cfaN] = sensorDetermineCFA(sensor);
%   [cfa,cfaN,mp] = sensorDetermineCFA; 
%   figure; image(cfaN); colormap(gray);
%
% See also: sensorColorOrder, sensorImageColorArray

% Programming notes
%
% In the future we could derive a color from the spectra using
% simple methods that we apply in scripts from time-to-time.
%
% --
% History - Until May 2010, the CFAnumbers would contain the value 4 at all
% positions with a cyan (c) filter.   The filters in every sensor must have
% distinct first letters.
%
% Until May 2010, the CFAnumbers referred to the position within this list,
% not the position of the sensor in the colorFilters columns. Because the
% field sensor.cfa.pattern uses numbers that refer to the column position
% in the colorFilters, and because we would like that field to work well
% with multispectral (> 3) imagers, we started to implement this change.
%
% Programming notes:
% Until May 2010, The color filter names and CFAletters mainly serve as a
% hint to the color appearance. The CFAnumbers numbers did NOT refer to the
% column in the sensor.color.filterSpectra, despite what was written here.
%
% The code used both letters and numbers to define each of these filters
% when plotting curves and generating and sensor image representations. 
%  
% Confusingly, there are really 3 different roles that in a better world
% would  correspond to filterNames, filterSpectra, filterColorHint. We use
% filterNames to serve both as names and as filterColorHints also, and
% this makes things more confusing than they should be. 
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('sensor'),
    sensor = vcGetObject('sensor');
    if isempty(sensor), error('no sensor defined'); end
end

rows = sensorGet(sensor,'rows');
cols = sensorGet(sensor,'cols');

pattern = sensorGet(sensor,'pattern');
if size(pattern,1)  ~= rows || size(pattern,2) ~= cols   
    % We have a cfa pattern that must match the sensor size
    blockRows = sensorGet(sensor,'unit block rows'); % size(pattern,1);
    blockCols = sensorGet(sensor,'unit block cols'); % size(pattern,2);
    
    cFactor = cols/blockCols;
    rFactor = rows/blockRows;
    if (round(cFactor) ~= cFactor) || (round(rFactor) ~= rFactor),
        error('Array size must be integer multiples of the size of the unit CFA block');
    end
else
    rFactor = 1; cFactor = 1;
end

CFAnumbers = repmat(pattern,rFactor,cFactor);

% Create the list of characters that are hints to the color appearance we
% should assign to each color filter. These hints can be useful for Matlab
% plotting routines and some ISET display routines. Get the letters from
% the first character of the filter name. 
filterColorLetters = sensorGet(sensor,'filterColorLetters'); 

% Figure out the filters that are OK with the first letter naming
% convention.
patternColors = sensorGet(sensor,'patternColors');

% Repeat the small block of letters to be the same size as the whole.
CFAletters = repmat(patternColors,rFactor,cFactor);

% Return the map
if nargout > 2
    [knownColorLetters, knownMap] = sensorColorOrder('string');
    nLetters = length(filterColorLetters);
    mp = zeros(nLetters,3);
    % Get this map from the known map
    for ii=1:nLetters
        idx = filterColorLetters(ii) == knownColorLetters;
        mp(ii,:) = knownMap(idx,:);
    end
end

end

