function obj = bipolarCompute(obj, os)
% Computes the responses of the bipolar object. 
% 
% The response is found by convolving the Gaussian spatial receptive field
% of the bipolar with each input image frame, convolving that signal with
% the temporal impulse response of the bipolar cell and thresholding.
% 
% 5/2016 JRG (c) isetbio team

% Spatial convolution
spatialResponseCenter = ieSpaceTimeFilter(os.coneCurrentSignal, obj.sRFcenter);
spatialResponseSurround = ieSpaceTimeFilter(os.coneCurrentSignal, obj.sRFsurround);

% Subsample to pull out individual bipolars
strideSubsample = size(obj.sRFcenter,1);
spatialSubsampleCenter = ieImageSubsample(spatialResponseCenter, strideSubsample);
spatialSubsampleSurround = ieImageSubsample(spatialResponseSurround, strideSubsample);

% Reshape for temporal convolution
szSubSample = size(spatialSubsampleCenter);
spatialSubsampleCenterRS = reshape(spatialSubsampleCenter,szSubSample(1)*szSubSample(2),szSubSample(3));
spatialSubsampleSurroundRS = reshape(spatialSubsampleSurround,szSubSample(1)*szSubSample(2),szSubSample(3));

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

obj.responseCenter = abs(bipolarOutputLinearCenter);
obj.responseSurround = zeros(size(bipolarOutputLinearSurround));

% % No - nonlinearity occurs in RGC computation after dot product with RGC RF
% % % Threshold
% % if ~isempty(obj.threshold)
% %     eZero = obj.threshold;
% %     obj.response = ieHwrect(bipolarOutputLinear,eZero);
% % else
% %     obj.response = bipolarOutputLinear;
% % end