function obj = osCompute(obj, sensor, varargin)
% Compute the linear filter response of the outer segments. 
%
%    obj = osCompute(obj, sensor, varargin)
%
% This converts isomerizations (R*) to outer segment current (pA). If the
% noiseFlag is set to true (1), this method adds noise to the current
% output signal. See Angueyra and Rieke (2013, Nature Neuroscience) for
% details.
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

if size(varargin) ~= 0
    if isfield(varargin{1,1},'convolutionType')
        cType = varargin{1,1}.convolutionType; % need to make this an input parameter!
    else
        cType = 0;
    end
else
    cType = 0;
end

% Remake filters incorporating the sensor to make them the correct sampling
% rate.
obj.matchSensor(sensor);

obj.patchSize = sensorGet(sensor,'width','meters'); % Patch of cone size

obj.timeStep  = sensorGet(sensor,'time interval','sec'); % Temporal sampling

% Find coordinates of L, M and S cones, get voltage signals.
cone_mosaic = sensorGet(sensor,'cone type');

% When we just use the number of isomerizations, this is consistent with
% the old coneAdapt function and validates.  
isomerizations = sensorGet(sensor,'photon rate');

% Get number of time steps.
% nSteps = sensorGet(sensor,'n time frames');
nSteps = size(sensor.data.volts,3);

% The next step is to convolve the 1D filters with the 1D isomerization
% data at each point in the cone mosaic. This code was adapted from the
% osLinearCone.m file by FR and NC.

initialState = osInit;
initialState.timeInterval = sensorGet(sensor, 'time interval');
initialState.Compress = 0; % ALLOW ADJUST - FIX THIS

% Place limits on the maxCur and prescribe the meanCur.

% See Angueyra and Rieke (2013, Nature Neuroscience)
maxCur = initialState.k * initialState.gdark^initialState.h/2;
meanCur = maxCur * (1 - 1 / (1 + 45000 / mean(isomerizations(:))));

% adaptedDataRS = osConvolve(obj, sensor, isomerizations, varargin);

[sz1, sz2, sz3] = size(isomerizations);
isomerizationsRS = reshape(isomerizations(:,:,1:sz3),sz1*sz2,nSteps);

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
    
    if ~isempty(cone_locations)
    isomerizationsSingleType = isomerizationsRS(cone_locations,:);
    
    % pre-allocate memory
    adaptedDataSingleType = zeros(size(isomerizationsSingleType));
    
    for y = 1:size(isomerizationsSingleType, 1)
        
        if ~cType     % Do full convolution with onset and offset transients
            tempData = conv(isomerizationsSingleType(y, :), FilterConeType);
        
        else % Do min length circular convolution in Fourier domain for steady state response only
            isoZP = [isomerizationsSingleType(y, :) zeros(1,-size(isomerizationsSingleType(y, :),2) + size(FilterConeType,1))];
            FilterConeTypeZP = [FilterConeType' zeros(1,size(isomerizationsSingleType(y, :),2)- size(FilterConeType,1))];
            tempData = ifft(fft(isoZP).*fft(FilterConeTypeZP));
        end 
        
        %  tempData = real(ifft(conj(fft(squeeze(isomerizationsSpec(x, y, :))) .* FilterFFT)));
        %  WHAT IS THE COMPRESS ABOUT?  Let's ASK NC
        if (initialState.Compress)
            tempData = tempData / maxCur;
            tempData = meanCur * (tempData ./ (1 + 1 ./ tempData)-1);
        else
            tempData = tempData - meanCur;
        end
        % NEED TO CHECK IF THESE ARE THE RIGHT INDICES

        if ~cType
            adaptedDataSingleType(y, :) = tempData([2:1+nSteps]);
        else
            
            adaptedDataSingleType(y, 1:length(tempData)-1) = tempData([2:end]);
            adaptedDataSingleType(y, length(tempData):end) = tempData(end);
        end
    end    
    % Check this size(isoRS,2) is right for circ conv
    adaptedDataRS(cone_locations,:) = adaptedDataSingleType(:,1:size(isomerizationsRS,2));  
    end
end

% % Reshape the output signal matrix.
adaptedData = reshape(adaptedDataRS,[sz1,sz2,sz3]);

% obj.coneCurrentSignal = adaptedData;
obj = osSet(obj, 'cone current signal', adaptedData);

% Add noise
% The osAddNoise function expects and input to be isomerization rate.
% This is handled properly because the params has the time sampling
% rate included.
if osGet(obj,'noiseFlag') == 1
    params.sampTime = sensorGet(sensor, 'time interval');
    ConeSignalPlusNoiseRS = osAddNoise(adaptedDataRS, params); 
    coneCurrentSignalPlusNoise = reshape(ConeSignalPlusNoiseRS,[sz1,sz2,nSteps]);
    % obj.coneCurrentSignalPlusNoise = reshape(ConeSignalPlusNoiseRS,[sz1,sz2,nSteps]);
    
    obj = osSet(obj, 'coneCurrentSignal', coneCurrentSignalPlusNoise);
end

end


