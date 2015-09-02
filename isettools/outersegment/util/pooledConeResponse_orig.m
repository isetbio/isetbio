function [pooledData] = pooledConeResponse(obj, sensor)

% pooledConeResponse: a utility function that computes the pooled response
% across the cone mosaic as a step in the ideal observer computation found
% in Hass, Horwitz, Angueyra, Lindbloom-Brown & Rieke, "Chromatic Detection
% from cone photorectpors to V1 neurons to behavior in rhesus monkeys,"
% (2015).

% Find coordinates of L, M and S cones.
cone_mosaic = sensorGet(sensor,'cone type');
[sz1, sz2] = size(cone_mosaic);

% Get cone current signal for each cone in the mosaic over time.
coneCurrent = obj.ConeCurrentSignal;

% Get number of time steps.
nSteps = size(coneCurrent, 3);

% The next step is to convolve the 1D filters with the 1D current
% data at each point in the cone mosaic. 

[sz1, sz2, sz3] = size(coneCurrent);
coneCurrentRS = reshape(coneCurrent(:,:,1:sz3),[sz1*sz2],nSteps);

totalIters = 400;
fprintf('\nGenerating pooled noisy responses:\n');
for iter = 1:totalIters
    
    fprintf('\b\b\b%02d%%', round(100*iter/totalIters));
    
for cone_type = 2:4
    % Pull out the appropriate 1D filter for the cone type.
    % Filter_cone_type = newIRFs(:,cone_type-1);
    switch cone_type
        case 2
            FilterConeType = obj.sConeFilter;
        case 3
            FilterConeType = obj.mConeFilter;
        case 4
            FilterConeType = obj.lConeFilter;
    end
    FilterConeType = (FilterConeType - mean(FilterConeType))./max(FilterConeType - mean(FilterConeType));
    
    
    % Only place the output signals corresponding to pixels in the mosaic
    % into the final output matrix.
    cone_locations = find(cone_mosaic==cone_type);
    
    coneCurrentRSnoisy = riekeAddNoise(coneCurrentRS);
    coneCurrentSingleType = (coneCurrentRSnoisy(cone_locations,:));
    
    if (ndims(coneCurrent) == 3)
        
        % pre-allocate memory
        adaptedDataSingleType = zeros(size(coneCurrentSingleType));
        
        for y = 1:size(coneCurrentSingleType, 1)
            noisySignal = squeeze((coneCurrentSingleType(y, :)));
            tempData = conv(noisySignal, FilterConeType);
            %        tempData = real(ifft(conj(fft(squeeze(coneCurrent(x, y, :))) .* FilterFFT)));

            adaptedDataSingleType(y, :) = tempData(1:nSteps);
             
        end
        
        %     elseif (ndims(coneCurrent) == 2)
        %
        %         % pre-allocate memory
        %         adaptedData = zeros(size(coneCurrent,1),timeBins);
        %
        %         for xy = 1:size(coneCurrent, 1)
        %             tempData = conv(squeeze(coneCurrent(xy, :)), Filter);
        %             if (initialState.Compress)
        %                 tempData = tempData / maxCur;
        %                 tempData = meanCur * (tempData ./ (1 + 1 ./ tempData)-1);
        %             else
        %                 tempData = tempData - meanCur;
        %             end
        %             adaptedData(xy, :) = tempData(1:timeBins);
        %         end
        %     end
        
        adaptedDataRS(cone_locations,:) = adaptedDataSingleType;
        pooledData(iter, cone_type-1) = mean(adaptedDataSingleType(:));
        
    end
    
    
end

% toc
end
