function obj = bipolarCompute(obj, os)
% Computes the responses of the bipolar object. 
% 
% The response is found by convolving the Gaussian spatial receptive field
% of the bipolar with each input image frame, convolving that signal with
% the temporal impulse response of the bipolar cell and thresholding.
% 
% 5/2016 JRG (c) isetbio team

% Spatial convolution
spatialResponse = ieSpaceTimeFilter(os.coneCurrentSignal, obj.sRF);

% Subsample to pull out individual bipolars
strideSubsample = size(obj.sRF,1);
spatialSubsample = ieImageSubsample(spatialResponse, strideSubsample);

% Reshape for temporal convolution
szSubSample = size(spatialSubsample);
spatialSubsampleRS = reshape(spatialSubsample,szSubSample(1)*szSubSample(2),szSubSample(3));

% Temporal convolution
bipolarOutputRS = convn(spatialSubsampleRS,obj.tIR','full');
% figure; plot(conv(spatialSubsampleRS(50,:),obj.tIR'))
% Back to original shape
bipolarOutputLinear = reshape(bipolarOutputRS,szSubSample(1),szSubSample(2),size(bipolarOutputRS,2));
% figure; plot(squeeze(bipolarOutputLinear(20,20,:)));
% Threshold
if ~isempty(obj.threshold)
    eZero = obj.threshold;
    obj.response = ieHwrect(bipolarOutputLinear,eZero);
else
    obj.response = bipolarOutputLinear;
end