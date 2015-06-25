function obj = outersegmentCompute(obj, sensor, param, varargin)
%Compute current output for outersegment class
% Implementation example here
% 
% Check bgvolts - make input
% 6/22/15 James Golden


% if notDefined('sensor'),  sensor = vcGetSelectedObject('sensor'); end
if notDefined('obj'), error('outersegment object must be passed'); end;

if ~exist('param','var') || isempty(param)
    error('Parameter field required.');
end
% if ~exist('val','var'),   error('Value field required.'); end;


param = ieParamFormat(param);  % Lower case and remove spaces
switch lower(param)
    
    % NEEDS UPDATE TO FRED'S LATEST PARAMETERS!
    
    case {1, 'linearfilter'}        
        % generate linear temporal filters for L, M, S responses
        [newIRFs, ~, ~] = cone_linear_filter();
        
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
        
    
    case {2, 'nonlinear1'} % find better name, more specific
        % See rieke<TAB> for explanations.
        
        % sensor = sensorCreate('human');
        % sensor = sensorSet(sensor,'pixel voltage swing',0.05);
        % v = rand(32,32,200)*sensorGet(sensor,'pixel voltage swing');
        % sensor = sensorSet(sensor,'volts',v);
        % vcAddObject(sensor); sensorWindow;
        % [~,adaptedData] = coneAdapt(sensor, 4);
        % s = 255 / max(abs(adaptedData(:)));
        %                 vcNewGraphWin
        %                 for ii=1:size(adaptedData,3)
        %                     imagesc(abs(s*adaptedData(:,:,ii)));
        %                     pause(0.02);
        %                 end
        % 
        % Fix this mplay(adaptedData,'I');
        
        p = riekeInit;
        expTime = sensorGet(sensor,'exposure time');
        sz = sensorGet(sensor,'size');
        
        % absRate = sensorGet(sensor,'absorptions per second');        
        pRate = sensorGet(sensor, 'photon rate');
                
        % Compute background adaptation parameters
        bgVolts = 1; % need to make this an input parameter!
        bgR = bgVolts / (sensorGet(sensor,'conversion gain')*expTime);
        
        initialState = riekeAdaptSteadyState(bgR, p, sz);
        obj.ConeCurrentSignal  = riekeAdaptTemporal(pRate, initialState);
        
        if obj.noiseflag == 1
            obj.ConeCurrentSignalPlusNoise = riekeAddNoise(obj.ConeCurrentSignal);
        end
                       
    otherwise
        error('unknown adaptation type');
        
end