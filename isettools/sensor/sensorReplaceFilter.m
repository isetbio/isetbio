function isa = sensorReplaceFilter(isa,whichFilter,newFilterFile)
%
%   isa = sensorReplaceFilter(isa,[whichFilter],[newFilterFile])
%
% Replace a color filter in the isa color filter list.  The data
% are read from a file.
%
% Examples:
%    [val,isa] = vcGetSelectedObject('ISA');
%    isa = sensorReplaceFilter(isa);
%
% Copyright Imageval, LLC 2002

if notDefined('isa'), [~, isa] = vcGetSelectedObject('ISA'); end
if notDefined('newFilterFile'), newFilterFile = []; end

if notDefined('whichFilter'), 
    filterNames = sensorGet(isa,'filterNames');
    replaceName = ieReadString('Enter filter name to replace:',filterNames{1}); 
    if isempty(replaceName), return; end;
    whichFilter = strmatch(replaceName,filterNames);     
    if isempty(whichFilter), return; end;
end

wave = sensorGet(isa,'wave');

[data,newFilterNames] = ieReadColorFilter(wave,newFilterFile);
nCols = size(data,2);
if nCols > 1
    str=sprintf('Data contains %.0f columns.  Choose one.', nCols);
    whichColumn = ieReadNumber(str);
    if isempty(whichColumn) | whichColumn < 1 | whichColumn > nCols, return; end
else
    whichColumn = 1;
end

data = data(:,whichColumn);
filterSpectra = sensorGet(isa,'filterspectra');
filterSpectra(:,whichFilter) = data;
isa = sensorSet(isa,'filterspectra',filterSpectra);

% Get a new filter name
newFilterName = ieReadString('Enter filter name to replace:',newFilterNames{whichColumn});
isa = sensorSet(isa,'editFilterNames',whichFilter,newFilterName);
% filterNames = ieSetFilterName(whichFilter,newFilterName,isa);
% isa = sensorSet(isa,'filternames',filterNames);

% We don't alter the isa.cfa.pattern.  We just changed the filter
% spectrum and the filter name.  The pattern stays the same.

return;

