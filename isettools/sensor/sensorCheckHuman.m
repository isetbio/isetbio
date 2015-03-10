function bool = sensorCheckHuman(sensor)
%Determine if this is a human sensor model
%
%   bool = sensorCheckHuman(sensor)
%
% (c) Imageval Consulting, LLC 2012

if notDefined('sensor'), error('sensor required'); end

if strfind(sensorGet(sensor,'name'),'human'), bool = 1; 
elseif isfield(sensor,'human'), bool = 1;
else bool = 0;
end

end