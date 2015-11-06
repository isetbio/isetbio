function obj = rgcCompute(obj, outerSegment, varargin)
% rgcCompute: a method of @rgc that computes the spiking output of the
% rgc mosaic to an arbitrary stimulus. See also @rgc/rgcCompute.m.
% 
%           rgc1 = rgcCompute(rgc1, os);
% 
% The spikes are computed using the post-spike filter for a refractory period. 
%  These computations are carried out using code from Pillow, Shlens, 
% Paninski, Sher, Litke, Chichilnisky, Simoncelli, Nature, 2008, licensed 
% for modification, which can be found at 
% 
% http://pillowlab.princeton.edu/code_GLM.html
% 
% Inputs: outersegment.
% 
% Outputs: the rgc object with spiking, raster and PSTH responses.
% 
% Example:
% rgc1 = rgcCompute(rgc1, identityOS);
% 
% (c) isetbio
% 09/2015 JRG

% Call to the superclass method rgcCompute.m
% Linear and nonlinear responses calculated in @rgc/rgcCompute
obj = rgcCompute@rgc(obj, outerSegment, varargin{:});


fprintf('     \n');
fprintf('Spike Generation:\n');
tic;

for cellTypeInd = 1:length(obj.mosaic)
    % Call Pillow code for spike response over N trials        
    spikeResponse = computeSpikesPSF(obj.mosaic{cellTypeInd});    
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'spikeResponse', spikeResponse);
    
    % Compute raster and PSTH responses for plotting using Pillow code.
    [raster, psth] = computePSTH(obj.mosaic{cellTypeInd});
    
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'rasterResponse',raster);
    obj.mosaic{cellTypeInd} = mosaicSet(obj.mosaic{cellTypeInd},'psthResponse',psth);
        
    clear spikeResponse raster psth
end
close;
toc;




