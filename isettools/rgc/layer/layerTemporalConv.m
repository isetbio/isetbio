function [linTS, componentTS] = layerTemporalConv(absBlurred, layer)
%% Calculate the linear component of the time series, linTS, for one layer
%  
%   [linTS componentTS] = layerTemporalConv(absBlurred, layer)
%
% Apply temporal impulse response filters to the spatially blurred
% receptive field component images. 
%
% Inputs:
%   absBlurred:  Spatially blurred and downsampled cone absorptions
%                (r,c,t,component)
%
% Returns: 
%   componentTS: the separate RF components (e.g., center and surround)
%     convolved with their own impulse response function
%   linTS:       the sum of the component TS, i.e., the total cone-driven
%     inputs to the RGCs 
%
% absBlurred is a 4D matrix.  The first two dimension are rows and columns
% of the cone matrix.  The third dimension is time.  The 4th dimension
% indicates which component of the receptive field (center, surround) is
% represented.  These data are obtained by spatial convolution of the cone
% absorptions with the rf-component (center/surround).
%
% See also: layerSpatialConv, rgcComputeSpikes
%
% Example:
%
% (c) 2010 Stanford Synapse Team 

%% Input format checking
if notDefined('absBlurred')
    error('Input blurred absorptions required');
end
if notDefined('layer') || ~isequal(class(layer),'rgcLayer')
    error('layer has to belong to the rgcLayer class');
end

%% Retrieve some values:
% rgcP = layer.parent;
% nT = rgcP.get('nFrames'); % number of frames
[r,c,nT,nComponents] = size(absBlurred);
nRGC = r*c;

%% retrieve locations and temporal response

% Retrieve layer parameters
inTR     = layer.get('inTR');    % Impulse responses per component

% Make sure the center and surround impulse responses each sum to 1.
% This should really be implemented in the creation phase.
inTR = inTR * diag([1/abs(sum(inTR(:,1))),1/abs(sum(inTR(:,2)))]);

% The absBlurred is spatially, but not yet temporally blurred.  The
% stimulus temporal samples are shorter than we need because we need to
% allow the convolution to extend beyond the number of stimulus temporal
% samples.  We initialize the time series with enough space for all the
% RGCs and for the extended time period.
padLength = size(inTR,1)-1;        % Comment this - it has to do with convolution (convn) below
tsLength  = nT + padLength;        % nT along with extra for impulse resp.
linTS     = zeros(nRGC,tsLength);  % RGC linear input time series

%% Create the linear time series per (rgc component)

% Perform the temporal impulse response convolution per RF component
componentTS = cell(1,nComponents);
for cc = 1:nComponents  % For the center and surround
    tmp = absBlurred(:, :, :, cc);
    tmp = reshape(tmp,r*c,nT);
    
    % What is inTR?  Input temporal response.
    componentTS{cc} = convn(tmp,inTR(:,cc)','full');
    
    % Add convolved signal for the center or surround component to the
    % total linTS
    linTS = linTS + componentTS{cc}; %componentTS{cc}(:,1:tsLength);
end

end