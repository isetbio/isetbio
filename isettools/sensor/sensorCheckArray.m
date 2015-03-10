function sensorCheckArray(ISA)
% Visual check of the color filter array pattern
%
%     sensorCheckArray(ISA)
%
% The routine produces an image that shows the color filter array pattern.  
%
% Example:
%   sensorCheckArray(vcGetObject('ISA'));
%        
% Copyright ImagEval Consultants, LLC, 2003.

% This seems to work with the new sensorDetermineCFA 5/10/2010

% The config field contains (x,y,color) information for the entire array
% If the ISA array is big, then just plot a portion of it
cfa = sensorDetermineCFA(ISA);

if size(cfa,1) > 50,  cfa = cfa(1:50,1:50); end

sensorImageColorArray(cfa);

return;
