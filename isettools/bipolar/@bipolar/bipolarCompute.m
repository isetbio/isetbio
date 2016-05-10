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

% Apply differentiator
%     obj.temporalDifferentiator;    % differentiator function

% Temporal convolution
bipolarOutputCenterRS = convn(spatialSubsampleCenterRS,obj.tIR','full');
bipolarOutputSurroundRS = convn(spatialSubsampleSurroundRS,obj.tIR','full');
% figure; plot(conv(spatialSubsampleRS(50,:),obj.tIR'))

% Back to original shape
bipolarOutputLinearCenter = reshape(bipolarOutputCenterRS,szSubSample(1),szSubSample(2),size(bipolarOutputCenterRS,2));
bipolarOutputLinearSurround = reshape(bipolarOutputSurroundRS,szSubSample(1),szSubSample(2),size(bipolarOutputSurroundRS,2));
% figure; plot(squeeze(bipolarOutputLinear(20,20,:)));

obj.responseCenter = bipolarOutputLinearCenter;
obj.responseSurround = bipolarOutputLinearSurround;

% % No - nonlinearity occurs in RGC computation after dot product with RGC RF
% % % Threshold
% % if ~isempty(obj.threshold)
% %     eZero = obj.threshold;
% %     obj.response = ieHwrect(bipolarOutputLinear,eZero);
% % else
% %     obj.response = bipolarOutputLinear;
% % end