function obj = osCompute(obj, sensor)
% osCompute: this method of @osIdentity passes on the cone isomerizations
% (R*) without any temporal filtering. This subclass is intended to be used
% for stimulus-referred retinal ganglion cell models.
%
% Inputs: the osIdentity object and the sensor object.
% 
% Outputs: the osIdentity object, with the cone outer segment current set 
% to the cone isomerizations and optionally a noisy version of the
% isomerizations.
% 
% 8/2015 JRG

% Get isomerization array to convert to current (pA).
isomerizations = sensorGet(sensor, 'photons');

% For the osIdentity object, there is no temporal filtering on the
% isomerizations.

adaptedData = isomerizations;
obj.ConeCurrentSignal = adaptedData./max(adaptedData(:));

% Add noise if the flag is set.
if obj.noiseFlag == 1
    params.sampTime = sensorGet(sensor, 'time interval');
    ConeSignalPlusNoiseRS = riekeAddNoise(adaptedDataRS, params); close;
    obj.ConeCurrentSignalPlusNoise = reshape(ConeSignalPlusNoiseRS,[sz1,sz2,nSteps]);
    
%     if size(varargin) ~= 0
%         if isfield(varargin{1,1},'offset')
%             obj.ConeCurrentSignalPlusNoise = obj.ConeCurrentSignalPlusNoise - obj.ConeCurrentSignalPlusNoise(:, :, nSteps) - varargin{1,1}.offset;
%         end
%     end
    
end


