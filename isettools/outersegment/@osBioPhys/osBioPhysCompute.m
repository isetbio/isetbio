function obj = osBioPhysCompute(obj, sensor, varargin)
% osBioPhysCompute: a method of @osBioPhys that computes the output
% response of the L, M and S cone outer segments. This converts
% isomerizations (R*) to outer segment current (pA). The differential
% equation model by Rieke is applied here.
% If the noiseFlag  property of the osLinear object is set to 1, this 
% method will add noise to the current output signal.
% 
% http://isetbio.github.io/isetbio/cones/adaptation%20model%20-%20rieke.pdf
% and 
% https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
% 
% 8/2015 JRG NC DHB
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());

%     if size(varargin{1,1})==0
%         params = cell(1,1);
%     else
%          params = (varargin{1,1});
%     end
    
    p = riekeInit;
    expTime = sensorGet(sensor,'exposure time');
    sz = sensorGet(sensor,'size');
    
    % absRate = sensorGet(sensor,'absorptions per second');
    pRate = sensorGet(sensor, 'photon rate');
    
    % Compute background adaptation parameters
%     if isfield(params{1,1},'bgVolts')
%         bgVolts = params{1,1}.bgVolts; % need to make this an input parameter!
%     else
%         bgVolts = 0;
%     end
    bgVolts = sensorGet(sensor,'adaptation offset');
    bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);
    
    initialState = riekeAdaptSteadyState(bgR, p, sz);
    
    initialState.timeInterval = sensorGet(sensor, 'time interval');
    obj.ConeCurrentSignal  = riekeAdaptTemporal(pRate, initialState);
    
%     if isfield(params{1,1},'dc')
%         nSamples = length(obj.ConeCurrentSignal);
%         obj.ConeCurrentSignal = obj.ConeCurrentSignal - obj.ConeCurrentSignal(:, :, nSamples);
%     end
    
    if obj.noiseFlag == 1
%         coneSamplingRate = 825;
        sampTime = 1/sensorGet(sensor, 'time interval');
        obj.ConeCurrentSignalPlusNoise = riekeAddNoise(obj.ConeCurrentSignal*sampTime)./sampTime;

%         obj.ConeCurrentSignalPlusNoise = riekeAddNoise(obj.ConeCurrentSignal, paramsNoise);        
        close;
        
    end
end

