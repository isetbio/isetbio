function obj = osCompute(obj, sensor, varargin)
% osLinearCompute: a method of @osLinear that computes the linear filter
% response of the L, M and S cone outer segments. This converts
% isomerizations (R*) to outer segment current (pA). If the noiseFlag
% property of the osLinear object is set to 1, this method will add noise
% to the current output signal.
%
% adaptedOS = osLinearCompute(adaptedOS, sensor);
%
% Inputs: the osLinear object, the sensor object and an optional parameters
% field. params.offest determines the current offset.
% 
% Outputs: the osLinear object, with the cone outer segment current and
% optionally a noisy version of the cone outer segment current.
% 
% 8/2015 JRG NC DHB


% Remake filters incorporating the sensor to make them the 
% correct sampling rate.
obj.filterKernel(sensor);

% Find coordinates of L, M and S cones, get voltage signals.
cone_mosaic = sensorGet(sensor,'cone type');

% Get isomerization array to convert to current (pA).
isomerizations = sensorGet(sensor, 'photon rate');

% For the osIdentity object, there is no temporal filtering on the
% isomerizations.

adaptedData = isomerizations;
obj.ConeCurrentSignal = adaptedData;

% Add noise if the flag is set.
if obj.noiseFlag == 1
    params.sampTime = sensorGet(sensor, 'time interval');
    ConeSignalPlusNoiseRS = riekeAddNoise(adaptedDataRS, params); close;
    obj.ConeCurrentSignalPlusNoise = reshape(ConeSignalPlusNoiseRS,[sz1,sz2,nSteps]);
    
    if size(varargin) ~= 0
        if isfield(varargin{1,1},'offset')
            obj.ConeCurrentSignalPlusNoise = obj.ConeCurrentSignalPlusNoise - obj.ConeCurrentSignalPlusNoise(:, :, nSteps) - varargin{1,1}.offset;
        end
    end
    
end


