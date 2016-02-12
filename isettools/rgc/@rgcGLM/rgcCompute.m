function obj = rgcCompute(obj, outerSegment, varargin)
%% Generate the linear and nonlinear responses for an rgc object 
% Computes the linear and nonlinear output of the rgc mosaic to an
% arbitrary stimulus. These computations are carried out using code from
% Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky, Simoncelli, Nature,
% 2008, which can be found at
% 
% http://pillowlab.princeton.edu/code_GLM.html
% 
% Inputs: the outersegment object.
% 
% Outputs: the rgc object with linear and nonlinear responses.
% 
% Example:
% rgc1 = rgcCompute(rgc1, identityOS);
% 
% (c) isetbio
% 09/2015 JRG

%% The superclass rgcCompute carries out convolution of the linear STRF:
obj = rgcCompute@rgc(obj, outerSegment, varargin{:});

% fprintf('     \n');
% fprintf('Spike Generation:\n');
% tic;
% for cellTypeInd = 1:length(obj.mosaic)
%     % Call Pillow code to compute spiking outputs for N trials
%     spikeResponse = computeSpikesGLM(obj.mosaic{cellTypeInd,1});   
%     obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'spikeResponse', spikeResponse);
%         
%     % Call Pillow code to compute rasters and PSTHs
%     [raster, psth] = computePSTH(obj.mosaic{cellTypeInd,1});
%     obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'rasterResponse',raster);
%     obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'psthResponse',psth);
%         
%     clear spikeResponse raster psth
% end
% close;
% toc;




