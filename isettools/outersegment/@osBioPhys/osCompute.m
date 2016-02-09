function obj = osCompute(obj, sensor, varargin)
% osCompute: a method of @osBioPhys that computes the output
% response of the L, M and S cone outer segments. This converts
% isomerizations (R*) to outer segment current (pA). The difference
% equation model by Rieke is applied here.
% If the noiseFlag  property of the osLinear object is set to 1, this
% method will add noise to the current output signal.
%
% http://isetbio.github.io/isetbio/cones/adaptation%20model%20-%20rieke.pdf
% and
% https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% Inputs: the osBioPhys object, the sensor object and an optional parameters
% field. params.offest determines the current offset.
%
% Outputs: the osBioPhys object, with the cone outer segment current and
% optionally a noisy version of the cone outer segment current.
%
% adaptedOS = osCompute(adaptedOS, sensor);
%
% 8/2015 JRG NC DHB


obj.coneSpacing = sensorGet(sensor,'width','um'); % Cone width

obj.coneSampling  = sensorGet(sensor,'time interval','sec'); % Temporal sampling

p = osInit;
expTime = sensorGet(sensor,'exposure time');
sz = sensorGet(sensor,'size');

% absRate = sensorGet(sensor,'absorptions per second');
pRate = sensorGet(sensor, 'photon rate');
nSteps = size(pRate, 3);
% Compute background adaptation parameters

if size(varargin) ~= 0
    if isfield(varargin{1,1},'bgVolts')
        bgVolts = varargin{1,1}.bgVolts; % need to make this an input parameter!
    else
        bgVolts = 0;
    end
else
    bgVolts = 0;
end
bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);


initialState = osAdaptSteadyState(bgR, p, sz);

initialState.timeInterval = sensorGet(sensor, 'time interval');

% obj.coneCurrentSignal  = osAdaptTemporal(pRate, initialState);
coneCurrentSignal  = osAdaptTemporal(pRate, initialState);

obj = osSet(obj, 'coneCurrentSignal', coneCurrentSignal);

if osGet(obj,'noiseFlag') == 1
    % obj.coneCurrentSignalPlusNoise = osAddNoise(obj.coneCurrentSignal);
    coneCurrentSignalPlusNoise = osAddNoise(coneCurrentSignal);
    % obj.coneCurrentSignalPlusNoise = osAddNoise(obj.coneCurrentSignal, paramsNoise);
    
    obj = osSet(obj, 'coneCurrentSignal', coneCurrentSignalPlusNoise);
    
end
end

