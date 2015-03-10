function sformatted = stringFormat(s)
% Converts s to a standard string format for rgcParameters set/get calls
%
%    sformatted = stringFormat(s)
%
% The string is sent to lower case and spaces are removed.
%
% ex:
%     stringFormat('Spatial sampling Rate')

warning('%s: Deprecated. Use ieParamFormat instead', mfilename);
if ~ischar(s), error('s should be a string'); end

% Lower case
sformatted = lower(s);

% Remove spaces
sformatted = strrep(sformatted,' ','');

end


