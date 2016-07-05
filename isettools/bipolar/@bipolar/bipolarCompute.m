function obj = bipolarCompute(obj, os, varargin)
% Computes the responses of the bipolar object. 
% 
% The outersegment input contains frames of cone mosaic signal at a
% particular time step. The bipolar response is found by first convolving
% the center and surround Gaussian spatial receptive fields of the bipolar
% cell within each cone signal frame. Then, the resulting signal is put
% through the weighted temporal differentiator in order to result in an
% impulse response that approximates the IR of the RGC.
% 
% Particular options that could be employed are rezeroing of the signal at
% the end of the temporal computation as well as rectification on the
% output signal.
% 
% 5/2016 JRG (c) isetbio team

%% parse input parameters
p = inputParser;
p.addRequired('obj', @(x) isa(x, 'bipolar'));
p.addRequired('os', @(x) isa(x, 'outerSegment'));

% parse
p.parse(obj, os, varargin{:});

%% Spatial filtering and subsampling
% Convolve spatial RFs over whole image, subsample to get evenly spaced
% mosaic.

% Get zero mean cone current signal
osSig = RGB2XWFormat(os.coneCurrentSignal);
osSig = bsxfun(@minus, osSig, mean(osSig, 2));
osSigZM = reshape(osSig, size(os.coneCurrentSignal));

% Compute spatial center/surround response
sCenter = convn(osSigZM, obj.sRFcenter, 'same');
sSurround = convn(osSigZM, obj.sRFsurround, 'same');

% Subsample to pull out individual bipolars
stride = size(obj.sRFcenter, 1);
sCenter = sCenter(1:stride:size(sCenter, 1), 1:stride:size(sCenter, 2), :);
sSurround = sSurround(1:stride:size(sSurround, 1), ...
    1:stride:size(sSurround, 2), :);

%% Temporal filtering
% Apply the weighted differentiator to the output of the spatial response

% Reshape for temporal convolution
[sCenter, r, c] = RGB2XWFormat(sCenter);
sSurround = RGB2XWFormat(sSurround);

% Repeat input matrix for filtering
sCenter = [repmat(sCenter(:,1), 1, 1) sCenter];
sSurround = [repmat(sSurround(:,1), 1, 1) sSurround];    

% load filters
if obj.filterType == 1  % average filter from measurement data
    if strcmpi(obj.cellType, 'offDiffuse')
        data = load('bipolarFilt.mat', 'bipolarOFFP');
        bipolarFilt = -data.bipolarOFFP(:)';
    else
        data = load('bipolarFilt.mat', 'bipolarONP');
        bipolarFilt = -data.bipolarONP(:)';
    end
    
    % Interpolate filter from data to get correct sample rate
    originalTimeStep = 0.008;
    bpLength = length(bipolarFilt);
    bipolarFilt = interp1(originalTimeStep:originalTimeStep:originalTimeStep*bpLength,bipolarFilt,obj.timeStep:obj.timeStep:originalTimeStep*bpLength);
    bipolarFilt(isnan(bipolarFilt)) = 0;
    
elseif obj.filterType == 2
    load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/irGLM.mat');
    if strcmpi(obj.cellType, 'offDiffuse')
        bipolarFilt = irGLM;
    else
        bipolarFilt = -irGLM;
    end
    
elseif obj.filterType == 3  % each temporal filter from the dataset
    data = load('/Users/james/Documents/MATLAB/isetbio misc/bipolarTemporal/bipolarFilt_200_OFFP_2013_08_19_6_all_linear.mat');
    bipolarFilt = (data.bipolarFiltMat(obj.cellLocation,:)');
end

% temporal filtering
tCenter = conv2(sCenter, bipolarFilt, 'same');
tSurround = conv2(sSurround, bipolarFilt, 'same');

% zero centering the output
tCenter = reshape(bsxfun(@minus, tCenter, mean(tCenter, 2)), r, c, []);
tSurround = reshape(bsxfun(@minus,tSurround,mean(tSurround, 2)), r, c, []);

%% Rectification
obj.responseCenter = obj.rectificationCenter(tCenter);
obj.responseSurround = obj.rectificationSurround(tSurround);

end