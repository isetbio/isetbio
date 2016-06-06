function obj = bipolarCompute(obj, os)
% Computes the responses of the bipolar object. 
% 
% The (x,y,t) input consists of "frames" which are the cone mosaic
% signal at a particular time step. The bipolar response is found by first
% convolving the center and surround Gaussian spatial receptive fields of
% the bipolar cell with each cone signal frame. Then, that resulting signal
% is put through the weighted temporal differentiator in order to result
% in an impulse response that approximates the IR of the RGC.
% 
% Particular options that could be employed are rezeroing of the signal at
% the end of the temporal computation as well as rectification on the
% output signal.
% 
% 5/2016 JRG (c) isetbio team

%% Spatial response
% Convolve spatial RFs over whole image, subsample to get evenly spaced
% mosaic.

% Get zero mean cone current signal
osSigRS = reshape(os.coneCurrentSignal, size(os.coneCurrentSignal,1)*size(os.coneCurrentSignal,2),size(os.coneCurrentSignal,3));
osSigRSZM = osSigRS - repmat(mean(osSigRS,2),1,size(osSigRS,2));
osSigZM = reshape(osSigRSZM,size(os.coneCurrentSignal));

% Set symmetric range
% osSigRS = reshape(os.coneCurrentSignal, size(os.coneCurrentSignal,1)*size(os.coneCurrentSignal,2),size(os.coneCurrentSignal,3));
% osSigRSrange = max(osSigRS(:))+85;
% osSigRSZM = osSigRSrange *( osSigRS - repmat(mean(osSigRS,2),1,size(osSigRS,2)));
% osSigZM = reshape(osSigRSZM,size(os.coneCurrentSignal));


% % % Spatial averaging over RGC blocks carried out here
% % % Need to move this before photoreceptors or else it doesn't matter
% rfSize = 39*size(obj.sRFcenter,1);
% numBlocksR = ceil(size(osSigZM,1)/(size(obj.sRFcenter,1)*rfSize));
% numBlocksC = ceil(size(osSigZM,2)/(size(obj.sRFcenter,1)*rfSize));
% for frInd = 1:size(osSigZM,3)
%     for blockIndR = 1:numBlocksR
%         for blockIndC = 1:numBlocksC
%             rindS = (blockIndR-1)*rfSize+1; rindE = (blockIndR)*rfSize; rindA = rindS:rindE;
%             rindV = rindA(rindA<size(osSigZM,1));
%             cindS = (blockIndC-1)*rfSize+1; cindE = (blockIndC)*rfSize; cindA = cindS:cindE;
%             cindV = cindA(cindA<size(osSigZM,2));
%             osSigSpAvg(rindV,cindV,frInd) = mean(mean(osSigZM(rindV,cindV,frInd)));
%         end
%     end
% end
            


% osSigZM(abs(osSigZM)>100) = 0;

% Spatial convolution
% spatialResponseCenter = ieSpaceTimeFilter(os.coneCurrentSignal-os.coneCurrentSignal(1,1,end), obj.sRFcenter);
% spatialResponseSurround = ieSpaceTimeFilter(os.coneCurrentSignal-os.coneCurrentSignal(1,1,end), obj.sRFsurround);

spatialResponseCenter = ieSpaceTimeFilter(osSigZM, obj.sRFcenter);
spatialResponseSurround = ieSpaceTimeFilter(osSigZM, obj.sRFsurround);
% spatialResponseCenter = ieSpaceTimeFilter(osSigSpAvg, obj.sRFcenter);
% spatialResponseSurround = ieSpaceTimeFilter(osSigSpAvg, obj.sRFsurround);


% Subsample to pull out individual bipolars
strideSubsample = size(obj.sRFcenter,1);
spatialSubsampleCenter = ieImageSubsample(spatialResponseCenter, strideSubsample);
spatialSubsampleSurround = ieImageSubsample(spatialResponseSurround, strideSubsample);

% 
% spatialSubsampleCenter = spatialResponseCenter;
% spatialSubsampleSurround = spatialResponseSurround;
%% Temporal response
% Apply the weighted differentiator to the output of the spatial
% computation.

% Reshape for temporal convolution
szSubSample = size(spatialSubsampleCenter);
spatialSubsampleCenterRS = reshape(spatialSubsampleCenter,szSubSample(1)*szSubSample(2),szSubSample(3));
spatialSubsampleSurroundRS = reshape(spatialSubsampleSurround,szSubSample(1)*szSubSample(2),szSubSample(3));

%%%% CHANGE THIS WHEN RUNNING DIFF!
obj.temporalDelay = 0;
% Zero pad to allow for delay
spatialSubsampleCenterRS = [repmat(spatialSubsampleCenterRS(:,1),1,(1e-3/os.timeStep)*obj.temporalDelay + 1).*ones(size(spatialSubsampleCenterRS,1),(1e-3/os.timeStep)*obj.temporalDelay + 1) spatialSubsampleCenterRS];
spatialSubsampleSurroundRS = [repmat(spatialSubsampleSurroundRS(:,1),1,(1e-3/os.timeStep)*obj.temporalDelay + 1).*ones(size(spatialSubsampleSurroundRS,1),(1e-3/os.timeStep)*obj.temporalDelay + 1) spatialSubsampleSurroundRS];    

% Apply the differentiator function.
% x = spatialSubsampleCenterRS;%-spatialSubsampleCenterRS(1);
% coneSig = obj.temporalConeW*x(:,2+obj.temporalDelay:end);
% coneDiff = obj.temporalConeDiffW*diff(x(:,1+obj.temporalDelay:end),1,2);
% figure; plot(coneSig+coneDiff); 
% hold on; plot(coneSig)
% figure; plot(x(:,2+obj.temporalDelay:end));
% figure; plot(diff(x(:,1+obj.temporalDelay:end),1,2));
% plot(coneDiff);

% % Differentiator
% bipolarOutputCenterRSLong = obj.temporalDifferentiator(spatialSubsampleCenterRS);
% bipolarOutputSurroundRSLong = obj.temporalDifferentiator(spatialSubsampleSurroundRS);

% % Convolve
% % DO THE CIRCULAR CONV
load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/bipolarFilt200_GLM.mat');
% bipolarOutputCenterRSLong = convn(spatialSubsampleCenterRS',bipolarFilt)';
% bipolarOutputSurroundRSLong =  convn(spatialSubsampleSurroundRS',bipolarFilt)';

% % % % Full ZP
% bipolarOutputCenterRSLongZP = [spatialSubsampleCenterRS zeros([size(spatialSubsampleCenterRS,1) size(bipolarFilt,1)])];
% bipolarOutputSurroundRSLongZP = [spatialSubsampleSurroundRS zeros([size(spatialSubsampleSurroundRS,1) size(bipolarFilt,1)])];
% bipolarFiltZP = repmat([bipolarFilt' zeros([size(bipolarFilt,2) size(spatialSubsampleCenterRS,2)])],size(spatialSubsampleCenterRS,1) ,1);

% % % Min ZP

if size(spatialSubsampleCenterRS,2) > size(bipolarFilt,1)
    bipolarOutputCenterRSLongZP = [spatialSubsampleCenterRS];% zeros([size(spatialSubsampleCenterRS,1) size(bipolarFilt,1)])];
    bipolarOutputSurroundRSLongZP = [spatialSubsampleSurroundRS];% zeros([size(spatialSubsampleSurroundRS,1)-size(bipolarFilt,1)])];
    bipolarFiltZP = repmat([bipolarFilt; zeros([-size(bipolarFilt,1)+size(spatialSubsampleCenterRS,2)],1)]',size(spatialSubsampleCenterRS,1) ,1);
else

    bipolarOutputCenterRSLongZP = ([spatialSubsampleCenterRS repmat(zeros([size(bipolarFilt,1)-size(spatialSubsampleCenterRS,2)],1)',size(spatialSubsampleCenterRS,1),1)]);
    
    bipolarOutputSurroundRSLongZP = ([spatialSubsampleSurroundRS repmat(zeros([size(bipolarFilt,1)-size(spatialSubsampleSurroundRS,2)],1)',size(spatialSubsampleSurroundRS,1),1)]);
    bipolarFiltZP = repmat(bipolarFilt',size(spatialSubsampleSurroundRS,1),1);
    
end
% 
bipolarOutputCenterRSLong = ifft(fft(bipolarOutputCenterRSLongZP').*fft(bipolarFiltZP'))';
bipolarOutputSurroundRSLong = ifft(fft(bipolarOutputSurroundRSLongZP').*fft(bipolarFiltZP'))';
%                 
% 
% 
bipolarOutputCenterRS = bipolarOutputCenterRSLong;%(:,1:end-(1e-3/os.timeStep)*obj.temporalDelay);
bipolarOutputSurroundRS = bipolarOutputSurroundRSLong;%(:,1:end-(1e-3/os.timeStep)*obj.temporalDelay);
% 
% bipolarOutputCenterRS =   [repmat(bipolarOutputCenterRSLong(:,1),1,round((1e-3/os.timeStep)*obj.temporalDelay)).*ones([size(bipolarOutputCenterRSLong,1)   round((1e-3/os.timeStep)*obj.temporalDelay)]) bipolarOutputCenterRSLong(:,1:end-(1e-3/os.timeStep)*obj.temporalDelay)];
% bipolarOutputSurroundRS =   [repmat(bipolarOutputSurroundRSLong(:,1),1,round((1e-3/os.timeStep)*obj.temporalDelay)).*ones([size(bipolarOutputSurroundRSLong,1)   round((1e-3/os.timeStep)*obj.temporalDelay)]) bipolarOutputSurroundRSLong(:,1:end-(1e-3/os.timeStep)*obj.temporalDelay)];

% Rezero
bipolarOutputCenterRSRZ = ((bipolarOutputCenterRS-repmat(mean(bipolarOutputCenterRS,2),1,size(bipolarOutputCenterRS,2))));
bipolarOutputSurroundRSRZ = ((bipolarOutputSurroundRS-repmat(mean(bipolarOutputSurroundRS,2),1,size(bipolarOutputSurroundRS,2))));
% figure; plot(conv(spatialSubsampleRS(50,:),obj.tIR'))

% Back to original shape
bipolarOutputLinearCenter = reshape(bipolarOutputCenterRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputCenterRS,2));
bipolarOutputLinearSurround = reshape(bipolarOutputSurroundRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputSurroundRS,2));
% figure; plot(squeeze(bipolarOutputLinear(20,20,:)));


%% Attach output to object
% % Bipolar rectification 
% obj.responseCenter = os.coneCurrentSignal;
% obj.responseSurround = zeros(size(os.coneCurrentSignal));

obj.responseCenter = (bipolarOutputLinearCenter);
obj.responseSurround = zeros(size(bipolarOutputLinearSurround));

% bipolarOutputRectifiedCenter = bipolarOutputLinearCenter.*(bipolarOutputLinearCenter>0);
% bipolarOutputRectifiedSurround = zeros(size(-bipolarOutputLinearSurround.*(bipolarOutputLinearSurround<0)));
% 
% % bipolarOutputRectifiedCenter = 1*abs(bipolarOutputLinearCenter);
% % bipolarOutputRectifiedSurround = zeros(size(bipolarOutputLinearSurround));
% 
% obj.responseCenter = bipolarOutputRectifiedCenter;
% obj.responseSurround = bipolarOutputRectifiedSurround;
