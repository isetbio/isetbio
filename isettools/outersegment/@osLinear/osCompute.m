function obj = osCompute(obj, sensor, varargin)
% Compute the linear filter response of the outer segments. 
%
%    obj = osCompute(obj, sensor, varargin)
%
% This converts isomerizations (R*) to outer segment current (pA). If the
% noiseFlag is set to true (1), this method adds noise to the current
% output signal.
%
% Inputs: 
%  obj: osLinear object
%  sensor: struct
%     optional parameters field. params.offest determines the current
%     offset. 
% 
% Outputs: 
%  osLinear object which includes the cone outer segment current 
%  optionally a noisy version of the cone outer segment current (noiseFlag)
% 
% 8/2015 JRG NC DHB

% Remake filters incorporating the sensor to make them the 
% correct sampling rate.
obj.initialize(sensor);

% Find coordinates of L, M and S cones, get voltage signals.
cone_mosaic = sensorGet(sensor,'cone type');

% Get isomerization array to convert to current (pA).
isomerizations = sensorGet(sensor, 'photons');

% Get number of time steps.
nSteps = sensorGet(sensor,'n time frames');
% size(sensor.data.volts,3);

% The next step is to convolve the 1D filters with the 1D isomerization
% data at each point in the cone mosaic. This code was adapted from the
% osLinearCone.m file by FR and NC.

initialState = osInit;
initialState.timeInterval = sensorGet(sensor, 'time interval');
initialState.Compress = 0; % ALLOW ADJUST - FIX THIS

% Place limits on the maxCur and prescribe the meanCur.
% This formula is based on a model from XXXX Rieke with YYYY and so forth
% Perhaps we should not bury it here but maybe it deserves its own function
% more because its importance.
%
maxCur = initialState.k * initialState.gdark^initialState.h/2;
meanCur = maxCur * (1 - 1 / (1 + 45000 / mean(isomerizations(:))));

[sz1, sz2, sz3] = size(isomerizations);
isomerizationsRS = reshape(isomerizations(:,:,1:sz3),[sz1*sz2],nSteps);

adaptedDataRS = zeros(size(isomerizationsRS));

% Do convolutions by cone type.
for cone_type = 2:4  % Cone type 1 is black (i.e., a hole in mosaic)
    
    % Pull out the appropriate 1D filter for the cone type.
    % Filter_cone_type = newIRFs(:,cone_type-1);
    switch cone_type
        case 4
            FilterConeType = obj.sConeFilter;
        case 3
            FilterConeType = obj.mConeFilter;
        case 2
            FilterConeType = obj.lConeFilter;
    end
    
    % Only place the output signals corresponding to pixels in the mosaic
    % into the final output matrix.
    cone_locations = find(cone_mosaic==cone_type);
    
    isomerizationsSingleType = isomerizationsRS(cone_locations,:);
    
    % pre-allocate memory
    adaptedDataSingleType = zeros(size(isomerizationsSingleType));
    
    for y = 1:size(isomerizationsSingleType, 1)
        
        tempData = conv(isomerizationsSingleType(y, :), FilterConeType);
        %  tempData = real(ifft(conj(fft(squeeze(isomerizationsSpec(x, y, :))) .* FilterFFT)));
        %  WHAT IS THE COMPRESS ABOUT?  Let's ASK NC
        if (initialState.Compress)
            tempData = tempData / maxCur;
            tempData = meanCur * (tempData ./ (1 + 1 ./ tempData)-1);
        else
            tempData = tempData - meanCur;
        end
        % NEED TO CHECK IF THESE ARE THE RIGHT INDICES
        adaptedDataSingleType(y, :) = tempData([2:1+nSteps]);
        
    end    
    
    adaptedDataRS(cone_locations,:) = adaptedDataSingleType;  
    
end

% Reshape the output signal matrix.
adaptedData = reshape(adaptedDataRS,[sz1,sz2,sz3]);
obj.ConeCurrentSignal = adaptedData;

% The user may send in a DC offset in the varargin list.
% We check for it here.  Is this OK?  Is it part of the referenced model?
if size(varargin) ~= 0
    if isfield(varargin{1,1},'offset')
        obj.ConeCurrentSignal = obj.ConeCurrentSignal - obj.ConeCurrentSignal(:, :, nSteps) - varargin{1,1}.offset;
    end
end

% Add noise
if obj.noiseFlag == 1
    params.sampTime = sensorGet(sensor, 'time interval');
    ConeSignalPlusNoiseRS = osAddNoise(adaptedDataRS, params); 
    obj.ConeCurrentSignalPlusNoise = reshape(ConeSignalPlusNoiseRS,[sz1,sz2,nSteps]);
    
    if size(varargin) ~= 0
        if isfield(varargin{1,1},'offset')
            obj.ConeCurrentSignalPlusNoise = obj.ConeCurrentSignalPlusNoise - obj.ConeCurrentSignalPlusNoise(:, :, nSteps) - varargin{1,1}.offset;
        end
    end
    
end


