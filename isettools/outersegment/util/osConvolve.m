function adaptedDataRS = osConvolve(obj, sensor, varargin)
% Convolve the cone linear filters with the isomerization signals.
% 
%   adaptedData = osConvolve(obj, isomerizations, varargin)
% 
% This does the convolution step from osCompute in @osLinear. See Angueyra 
% and Rieke (2013, Nature Neuroscience) for details. 
% 
% Inputs: 
%  obj: osLinear object
%  sensor: struct containing isomerization signal
%  optional parameters field: the compress flag restricts the max value of
%       the current output.
% 
% Outputs: 
%   adaptedData: the cone current signal (pA).
% 
% 8/2015 JRG NC DHB

% Remake filters incorporating the sensor to make them the 
% correct sampling rate.
% obj.initialize(sensor);

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

% See Angueyra and Rieke (2013, Nature Neuroscience) for details
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
        adaptedDataSingleType(y, :) = tempData([2:1+sz3]);
        
    end    
    
    adaptedDataRS(cone_locations,:) = adaptedDataSingleType;  
    
end

% % Reshape the output signal matrix.
% adaptedData = reshape(adaptedDataRS,[sz1,sz2,sz3]);