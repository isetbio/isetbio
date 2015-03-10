function sformatted = ieParamFormat(s)
% Converts s to a standard ISET parameter format  (lower case, no spaces)
%
%    sformatted = ieParamFormat(s)
%
% The string is sent to lower case and spaces are removed.
%
% Example:
%     ieParamFormat('Exposure Time')
%
% See also: sensorGet, and so forth
%
% Copyright ImagEval Consultants, LLC, 2010

if ~ischar(s), error('s has to be a string'); end

% Lower case
sformatted = lower(s);

% Remove spaces
sformatted = strrep(sformatted,' ','');

end


