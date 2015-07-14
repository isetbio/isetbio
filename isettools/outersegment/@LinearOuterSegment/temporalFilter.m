function temporalFilter(obj, sensor, param, varargin)
    
    fprintf('<strong>\n%s:\n\t%s()\n</strong>', class(obj), mfilename());
    
    %
    % call the parent class method to compute the filter
    % obj.computeFilter();
    % 
    % compute a linear filtering operation
    % obj.outputSignal = conv(obj.filterKernel, obj.inputSignal);
    
    
    
    % generate linear temporal filters for L, M, S responses
    [newIRFs, ~, ~] = cone_linear_filter();
%     obj.sConeFilter = newIRFs(:,2);
%     obj.mConeFilter = newIRFs(:,3);
%     obj.lConeFilter = newIRFs(:,4);
    
    % find coordinates of l, m, s cones, get voltage signals
    cone_mosaic = sensorGet(sensor,'cone type');
    [sz1, sz2] = size(cone_mosaic);
    
    % get isomerization array to convert to current (pA)
    isomerizations = sensorGet(sensor, 'photons');
    
    % get number of time steps
    nSteps = size(sensor.data.volts,3);
    
    for cone_type = 2:4
        % create rows X cols X time matrix of temporal filters
        Filter_cone_type = newIRFs(:,cone_type-1);
        Filter_block = repmat(fft(Filter_cone_type(1:nSteps)'),[1 size(sensor.data.volts,1) size(sensor.data.volts,2)]);
        Filter_block2 = reshape(Filter_block,size(sensor.data.volts));
        
        % filter isomerizations matrix
        
        %  MAKE THIS GENERAL
        coneSamplingRate = 825; % samples per second
        obj.ConeCurrentSignal = real(ifft((Filter_block2) .* fft(isomerizations,[],3),[],3)) / coneSamplingRate;
        
        % reshape signal matrix
        cone_locations = find(cone_mosaic==cone_type);
        ConeSignal_rs = reshape(obj.ConeCurrentSignal,[sz1*sz2],nSteps);
        ConeSignalFinal_rs(cone_locations,:) = ConeSignal_rs(cone_locations,:);
        % obj.ConeSignalFinalCell{cone_type-1} = ConeSignal_rs(cone_locs,:);
        ConeCurrentSignalCell{cone_type-1} = ConeSignal_rs(cone_locations,:);
        
    end
    
    
    % add noise if flag is set
    if obj.noiseflag == 1
        
        % rescale by sampling rate, add noise
        % CHECK IF RESCALE IS CORRECT
        ConeSignalPlusNoise_rs = riekeAddNoise(ConeSignalFinal_rs*coneSamplingRate)./coneSamplingRate;
        
        obj.ConeCurrentSignalPlusNoise = reshape(ConeSignalPlusNoise_rs,[sz1,sz2,nSteps]);
        % reshape
        for cone_type=2:4
            
            cone_locations = find(cone_mosaic==cone_type);
            ConeCurrentSignalPlusNoiseCell{cone_type-1} = ConeSignalPlusNoise_rs(cone_locations,:);
        end
        
    end
    
end

