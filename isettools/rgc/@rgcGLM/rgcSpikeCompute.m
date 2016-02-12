function obj = rgcSpikeCompute(obj, varargin)
%% Generate spikes from rgcGLM object response
% Computes the spiking output of the rgc mosaic to an arbitrary stimulus.
% These computations are carried out using code from Pillow, Shlens,
% Paninski, Sher, Litke, Chichilnisky, Simoncelli, Nature, 2008, licensed
% for modification, which can be found at
% 
% http://pillowlab.princeton.edu/code_GLM.html
% 
% Inputs: the rgc object with linear and nonlinear responses computed using
%       rgcCompute.
% 
% Outputs: the rgc object with spiking responses.
% 
% Example:
% rgc1 = rgcSpikeCompute(rgc1);
% 
% (c) isetbio
% 09/2015 JRG

%% Call the Pillow code to generate spikes for the whole mosaic using the coupled GLM
% The superclass rgcCompute carries out convolution of the linear STRF:
% obj = rgcCompute@rgc(obj, outerSegment, varargin{:});

% fprintf('     \n');
% fprintf('Spike Generation:\n');
% tic;
for cellTypeInd = 1:length(obj.mosaic)
    % Call Pillow code to compute spiking outputs for N trials
    spikeResponse = computeSpikesGLM(obj.mosaic{cellTypeInd,1});   
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'spikeResponse', spikeResponse);
        
    % Call Pillow code to compute rasters and PSTHs
%     [raster, psth] = computePSTH(obj.mosaic{cellTypeInd,1});
%     obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'rasterResponse',raster);
%     obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'psthResponse',psth);
        
    clear spikeResponse raster psth
end
% close;
% toc;




