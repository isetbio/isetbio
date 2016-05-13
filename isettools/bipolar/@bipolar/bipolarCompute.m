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

% Spatial convolution
spatialResponseCenter = ieSpaceTimeFilter(os.coneCurrentSignal, obj.sRFcenter);
spatialResponseSurround = ieSpaceTimeFilter(os.coneCurrentSignal, obj.sRFsurround);

% Subsample to pull out individual bipolars
strideSubsample = size(obj.sRFcenter,1);
spatialSubsampleCenter = ieImageSubsample(spatialResponseCenter, strideSubsample);
spatialSubsampleSurround = ieImageSubsample(spatialResponseSurround, strideSubsample);

%% Temporal response
% Apply the weighted differentiator to the output of the spatial
% computation.

% Reshape for temporal convolution
szSubSample = size(spatialSubsampleCenter);
spatialSubsampleCenterRS = reshape(spatialSubsampleCenter,szSubSample(1)*szSubSample(2),szSubSample(3));
spatialSubsampleSurroundRS = reshape(spatialSubsampleSurround,szSubSample(1)*szSubSample(2),szSubSample(3));

% Apply the differentiator function.
bipolarOutputCenterRS = obj.temporalDifferentiator(spatialSubsampleCenterRS);
bipolarOutputSurroundRS = obj.temporalDifferentiator(spatialSubsampleSurroundRS);

% Rezero
bipolarOutputCenterRSRZ = ((bipolarOutputCenterRS-repmat(mean(bipolarOutputCenterRS,2),1,size(bipolarOutputCenterRS,2))));
bipolarOutputSurroundRSRZ = ((bipolarOutputSurroundRS-repmat(mean(bipolarOutputSurroundRS,2),1,size(bipolarOutputSurroundRS,2))));
% figure; plot(conv(spatialSubsampleRS(50,:),obj.tIR'))

% Back to original shape
bipolarOutputLinearCenter = reshape(bipolarOutputCenterRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputCenterRS,2));
bipolarOutputLinearSurround = reshape(bipolarOutputSurroundRSRZ,szSubSample(1),szSubSample(2),size(bipolarOutputSurroundRS,2));
% figure; plot(squeeze(bipolarOutputLinear(20,20,:)));

%% Attach output to object
% Bipolar rectification 
obj.responseCenter = (bipolarOutputLinearCenter);
obj.responseSurround = zeros(size(bipolarOutputLinearSurround));


