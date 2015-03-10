function ISA = sensorClearData(ISA)
%Clear data and noise fields stored in the sensor array.
%
%   ISA = sensorClearData(ISA)
%
% When parameters change and data are no longer consistent, we clear the
% data and various stored noise image fields.  
%
% Copyright ImagEval Consultants, LLC, 2003.

if checkfields(ISA,'data'),           ISA = sensorSet(ISA,'data',[]); end
if checkfields(ISA,'offsetFPNimage'), ISA = sensorSet(ISA,'offsetFPNimage',[]); end
if checkfields(ISA,'gainFPNimage'),   ISA = sensorSet(ISA,'gainFPNimage',[]); end
if checkfields(ISA,'colOffset'),      ISA = sensorSet(ISA,'coloffsetfpnvector',[]); end
if checkfields(ISA,'colGain'),        ISA = sensorSet(ISA,'colgainfpnvector',[]); end

return;