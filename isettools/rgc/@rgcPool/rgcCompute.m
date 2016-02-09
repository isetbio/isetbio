function obj = rgcCompute(obj, outerSegment, varargin)
% rgcCompute: a method of @rgc that computes the spiking output of the
% rgc mosaic to an arbitrary stimulus.
%       
%           rgc1 = rgcCompute(rgc1, os);
% 
% The responses for each mosaic are computed one at a time. For a given
% mosaic, first the spatial convolution of the center and surround RFs are
% calculated for each RGB channel, followed by the temporal responses for
% the center and surround and each RGB channel. This results in the linear
% response.
% 
% Next, the linear response is put through the generator function. The
% nonlinear response is the input to a function that computes spikes with
% Poisson statistics. For the rgcGLM object, the spikes are computed using
% the recursive influence of the post-spike and coupling filters between
% the nonlinear responses of other cells. These computations are carried
% out using code from Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky,
% Simoncelli, Nature, 2008, licensed for modification, which can be found
% at 
% 
% http://pillowlab.princeton.edu/code_GLM.html
% 
% See @rgcGLM/rgcCompute.m for the implementation.
% 
% Outline:
% 1. Normalize stimulus
% 2. Compute linear response
%       - spatial convolution
%       - temporal convolution
% 3. Compute nonlinear response
% [spiking responses are calculated by subclass versions of rgcCompute]
% 
% Inputs: rgc object, outersegment object.
% 
% Outputs: the rgc object with responses.
% 
% Example:
% rgc1 = rgcCompute(rgc1, identityOS);
% 
% (c) isetbio
% 09/2015 JRG
%%

narginchk(0, Inf);
p = inputParser; p.CaseSensitive = false; p.FunctionName = mfilename;

% addRequired(p,'rgc',@isa();
addParameter(p,'mosaicType','',@ischar);

p.parse(varargin{:});
% 
% [obj.row, obj.col] = size(osGet(outerSegment,'coneCurrentSignal'));
% obj.spacing = osGet(outerSegment,'coneSpacing'); % Cone width
% obj.timing  = osGet(outerSegment,'coneSampling'); % Temporal sampling

if isempty(obj.mosaic{1})
    for cellTypeInd = 1:5%length(obj.mosaic)
        obj = rgcMosaicCreate(obj, 'mosaicType', p.Results.mosaicType);
    end
end

%%

% Get stimulus from outer segment ojbect
% Set zero mean and normalize max value of stimulus
if isa(outerSegment,'osIdentity')
    spTempStim = osGet(outerSegment, 'rgbData');
%     spTempStim = spTempStim - mean(spTempStim(:));% - 0.5;%1/sqrt(2);
%     spTempStim = 2*spTempStim./max(abs(spTempStim(:)));    
%     spTempStim = spTempStim - mean(spTempStim(:));% - 0.5;%1/sqrt(2);

    range = max(spTempStim(:))-min(spTempStim(:));
    
    if isa(obj,'rgcPhys'); 
        spTempStim = spTempStim./range; 
    else
        spTempStim = 10*spTempStim./range;
    end;
    
elseif isa(outerSegment,'osLinear')||isa(outerSegment,'osBioPhys')
    % This is after temporal processing - correct to set zero mean?
    spTempStim = osGet(outerSegment, 'coneCurrentSignal');    
    spTempStim = spTempStim - mean(spTempStim(:));
    spTempStim = 5*spTempStim./max(abs(spTempStim(:)));
end

nSamples = size(spTempStim,3);

for cellTypeInd = 1:length(obj.mosaic)
    
    rfSize = size(obj.mosaic{cellTypeInd}.sRFcenter{1,1});
    
    % Given separable STRF, convolve 2D spatial RF with each frame     
    [spResponseCenter, spResponseSurround] = spConvolve(obj.mosaic{cellTypeInd,1}, spTempStim);
    spResponse = (cellfun(@minus,spResponseCenter,spResponseSurround,'uniformoutput',false));
    spResponseVector = cellfun(@(x) squeeze(mean(mean(x,2),1)),spResponse,'uniformoutput',false);
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'linearResponse', spResponseVector);
    
%     if isa(outerSegment,'osIdentity')
%         % Then convolve output of spatial convolution with the temporal impulse response
% %         spResponseCenter = spTempStim; spResponseSurround = zeros(size(spTempStim));
%         [fullResponse, nlResponse] = fullConvolve(obj.mosaic{cellTypeInd,1}, spResponseCenter, spResponseSurround);
%     
%     elseif isa(outerSegment,'osLinear')||isa(outerSegment,'osBioPhys')
%         % Unless the os object has already applied temporal processing,
%         % then take output of spatial convolution as full output.
%         % Take difference between center and surround outputs of spatial
%         % convolution:
%         strfResponse = cellfun(@minus, spResponseCenter, spResponseSurround,'un',0);
%         % Find the mean over the strf response for each temporal frame:
%         strfResponseRS = cellfun(@(x) reshape(x,rfSize(1)*rfSize(2),nSamples),strfResponse,'un',0);
%         fullResponse = cellfun(@mean, strfResponseRS,'un',0);
%         if ~isa(obj, 'rgcLinear'); nlResponse = cellfun(obj.mosaic{cellTypeInd}.generatorFunction,fullResponse,'un',0); end;
%     end
%        
%     obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'linearResponse', fullResponse);
%         
%     % Set the nonlinear response for every rgc subclass except rgcLinear
%     if ~isa(obj, 'rgcLinear');
%         obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'nlResponse', nlResponse);
%     end
            
    clear spResponseCenter spResponseSurround fullResponse nlResponse 
end