function rgcComputeSpikes(rgcP)
% Converts cone volts to a spike time series
%
% The parameters are in the RGC parameter (rgcP) class.
%
% The input voltages from the cone mosaic are spatially and temporally
% convolved by the RGC spatio-temporal receptive fields to compute the
% linear time series (linTS).
%
% Then spikes are computed using a noisy process that includes feedback of
% cells onto themselves as well as network connection signals between the
% cells.
%
% Inputs :
%   rgcP:    rgcParameters, default is without any layer, to create
%            parameters with one 'ON' layer,
%
% Outputs :
%  The computed spikes, rgc time series, and linear time series are each
%  part of the layer structures in the rgcP.
% 
%  We get access to the computed values from each layer by calling
%  layer.get('linear time series'), or layer.get('spikes').  
%  For examples, see rgcVisualize.
%
% Example:
%   See the example of how to call this in s_rgcCones2RGC.m
%
% See also:  rgcVisualize
%
% (c) 2009 Stanford VISTA Team

%% Check inputs
if notDefined('rgcP') || ~isequal(class(rgcP),'rgcParameters')
    error('rgcP of rgcParameters class is required.');
end

nL = rgcP.get('nLayers');
if nL==0, error('There are no layers, use rgcP.addLayer() to add one'); end

%% Assigning and preprocessing the data

% You might check that these are OK
% rgcVisualize('cone voltages')

% Add a few frames of noise in the data to initialize the network in a
% random fashion (~steady state). If we use precomputed linTS, then don't
% add noise to the absorptions, it messes stuff up. Since the precomputed
% linTS should have been computed with the same number of noise frames, we
% will still have nFrames+nNoise spikes, and they will be trimmed by
% rgcRemoveNoise, called at the end of rgcComputeSpikesNoiseInit. - EC

% rgcP = rgcAddNoise(rgcP);

%% Computation of the layers

% lot of if else to optimise the code depending on what you want to do, but
% it is almost the same thing everytime.
% nT = rgcP.get('nFrames');
% rgcTS       = cell(1,nL);   % Final time series of cells in each layer
% spikes      = cell(1,nL);

for ii = 1:nL  % For each layer
    
    thisLayer = rgcP.getLayer(ii); % Get the parameters for this layer
    cellType = thisLayer.get('name');
    fprintf('[%s]: Layer %d (%s)\n',mfilename,ii,cellType);
    
    switch stringFormat(cellType)
        case {'onmidget','offmidget'}
            [spikes, rgcTS] = rgcComputeMidgets(thisLayer);
        case {'onparasol','offparasol'}
            [spikes, rgcTS] = rgcComputeParasol(thisLayer);
        case {'smallbistratified'}
            [spikes, rgcTS] = rgcComputeSmallBistratified(thisLayer);

        otherwise
            error('Unknown cell type %s\n',cellType)
    end
    
    
    % Reshape the spike array. The reshaping accounts for the grid size and
    % number of temporal samples. Then the permuting puts it into the
    % format that other parts of the rgc code expects for rows and columns
    % in the image.
    % nT = size(rgcTS,2);
    % spikes{ii}   = uint8(reshape(spikes{ii}, gS(1), gS(2), nT));
    
    % Clip the linear time series, the RGC time series, and the spikes to
    % match the duration of cone time series.
    coneFrames = rgcGet(rgcP,'nframes');
    spikes = spikes(:,1:coneFrames);
    %     gS = thisLayer.get('gridSize');
    %     spikes = uint8(reshape(spikes, gS(1), gS(2), coneFrames));

    linTS = thisLayer.get('linear time series');
    linTS = linTS(:,1:coneFrames);
    rgcTS = rgcTS(:,1:coneFrames);
    
    thisLayer.set('linear time series', linTS);
    thisLayer.set('rgc voltage time series',rgcTS);
    thisLayer.set('spike time series',spikes);
    
end

end


%% Computational methods for different RGC types from cone absorptions

%% Midget system
function [spikes, rgcTS] = rgcComputeMidgets(thisLayer)
%
% This is a linear followed by threshold method.
% We will implement different models here.
%
% Models we will implement:
%    Linear-nonlinear, threshold
%    Rectification ideas.
%    Coupling experiments and ideas
%
% Accumulate References here
%

connections = layersConnectionMatrix(thisLayer);
thisLayer.set('Connection Matrix',connections);

% Moment by moment spatial blurring of the absorptions
absBlurred = layerSpatialConv(thisLayer);
% absBlurred is 4D: (x,y,timeFrame,[center or surr component])
% The entries of absBlurred are in volts.

% Have a look at the spatially blurred data.
% vcNewGraphWin;
% imagesc(absBlurred(:,:,round(size(absBlurred,3)/2),1)); 
% mean(mean(absBlurred(:,:,round(size(absBlurred,3)/2),1)))

% From the spatially blurred absorptions, add the temporal impulse response
% for both center and surround
[linTS, componentTS] = layerTemporalConv(absBlurred, thisLayer);
thisLayer.set('lints',linTS);
thisLayer.set('RFcomponentTS', componentTS);

% rgcVisualize('Linear Time Series',rgcP);
% rgcVisualize('Linear Time Series image',rgcP);

% After we get the time series back, we might apply an 'adaptation'
% function to keep the signals in proper RGC voltage range.  This is
% analogous to performing cone adaptation - like contrast adaptation.
% We should also truncate the time series by vSwing, right?
%
linTS = layerAdapt(thisLayer);
% Maybe this should not over write the linTS, but be adaptTS.
thisLayer.set('lints',linTS);

% From the full linear time series (blurred in space and time) make spiking
% decisions
[spikes rgcTS] = layerComputeSpikes(thisLayer);



end

%% Parasol system
function [spikes, rgcTS] = rgcComputeParasol(thisLayer)
%
% This is just a copy of the midget, and for now we run this like the
% midget but with different parameters.  In the future, we will write a
% real model.
%
% This should include Shapley nonlinear subunit stuff, for example.
%

connections = layersConnectionMatrix(thisLayer);
thisLayer.set('Connection Matrix',connections);

absBlurred = layerSpatialConv(thisLayer);
% absBlurred is 4D: (x,y,timeFrame,[center or surr component])
% The entries of absBlurred are in volts.

% Have a look at the spatially blurred data.
% vcNewGraphWin;
% imagesc(absBlurred(:,:,round(size(absBlurred,3)/2),1));  % Layer 1
% mean(mean(absBlurred(:,:,round(size(absBlurred,3)/2),1)))
% imagesc(absBlurred(:,:,round(size(absBlurred,3)/2),2));  % Layer 2

[linTS, componentTS] = layerTemporalConv(absBlurred, thisLayer);
thisLayer.set('lints',linTS);
thisLayer.set('RFcomponentTS', componentTS);

% rgcVisualize('Linear Time Series',rgcP);
% rgcVisualize('Linear Time Series image',rgcP);

% After we get the time series back, we might apply an 'adaptation'
% function to keep the signals in proper RGC voltage range.  This is
% analogous to performing cone adaptation - like contrast adaptation.
% We should also truncate the time series by vSwing, right?
%
linTS = layerAdapt(thisLayer);
% Maybe this should not over write the linTS, but be adaptTS.
thisLayer.set('lints',linTS);

% Finally, use the linear time series and calculate the spikes
[spikes rgcTS] = layerComputeSpikes(thisLayer);

% Reshape the spikes. The reshaping accounts for the grid size and
% number of temporal samples. Then the permuting puts it into the
% format that other parts of the rgc code expects for rows and columns
% in the image.
% nT = size(rgcTS,2);
% gS = thisLayer.get('gridSize');
% 
% spikes   = uint8(reshape(spikes, gS(1), gS(2), nT));
% thisLayer.set('rgc voltage time series',rgcTS);
% thisLayer.set('spike time series',spikes);

end

%% Small bistratified
function [spikes, rgcTS] = rgcComputeSmallBistratified(thisLayer)
%
%  For the moment, it is the standard linear/threshold.  All the noise  We
%  is still in the cone absorptions here.
%

connections = layersConnectionMatrix(thisLayer);
thisLayer.set('Connection Matrix',connections);

absBlurred = layerSpatialConv(thisLayer);
% absBlurred is 4D: (x,y,timeFrame,[center or surr component])
% The entries of absBlurred are in volts.

% Have a look at the spatially blurred data.
% vcNewGraphWin;
% imagesc(absBlurred(:,:,round(size(absBlurred,3)/2),1));  % Layer 1
% mean(mean(absBlurred(:,:,round(size(absBlurred,3)/2),1)))
% imagesc(absBlurred(:,:,round(size(absBlurred,3)/2),2));  % Layer 2

[linTS,   componentTS] = layerTemporalConv(absBlurred, thisLayer);
thisLayer.set('lints',linTS);
thisLayer.set('RFcomponentTS', componentTS);

% rgcVisualize('Linear Time Series',rgcP);
% rgcVisualize('Linear Time Series image',rgcP);

% After we get the time series back, we might apply an 'adaptation'
% function to keep the signals in proper RGC voltage range.  This is
% analogous to performing cone adaptation - like contrast adaptation.
% We should also truncate the time series by vSwing, right?
%
linTS = layerAdapt(thisLayer);
% Maybe this should not over write the linTS, but be adaptTS.
thisLayer.set('lints',linTS);

% Finally, use the linear time series and calculate the spikes
[spikes rgcTS] = layerComputeSpikes(thisLayer);

% Reshape the spikes. The reshaping accounts for the grid size and
% number of temporal samples. Then the permuting puts it into the
% format that other parts of the rgc code expects for rows and columns
% in the image.
% gS = thisLayer.get('gridSize');
% nT = size(rgcTS,2);
% spikes   = uint8(reshape(spikes, gS(1), gS(2), nT));
% thisLayer.set('rgc voltage time series',rgcTS);
% thisLayer.set('spike time series',spikes);

end

