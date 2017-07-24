function [rgcL, nTrialsSpikes] = compute(rgcL, varargin)
% @RGCLAYER.COMPUTE - Computes the rgc mosaic responses to an input
%
%   @rgcLayer.compute(varargin)
%
% Computes the continuous (linear) and then spike responses for each of the
% mosaics within the inner retina layer object.
%
% Required inputs
%
% Optional inputs
%   'nTrialsSpikes' -  Multiple trials case
%
% Compute the responses for each mosaic in the layer. First, a space-time
% separable linear response is computed. This stage of the computation is
% stored in 'responseLinear'. This is managed in computeSeparable.  There
% is no noise added in the linear part.
%
% At present, the temporal response is set to an impulse because of the way
% the bipolar tIR is set.  We need to deal with this.  Also, the 'dt' of
% the RGC should be inherited form the dt of the bipolar calculation.  I
% don't see that anywhere (BW).  Should be in the
% irComputeLinearSTSeparable routine.
%
% The spikes are computed irComputeSpikes routine. The spiking can have a
% random element.  So, we may run the conversion from linear to spikes
% multiple times, effectively producing spike rasters.
%
% Outputs:
%  ir: the inner retina object with responses attached to each mosaic
%  nTrialsSpikes:  (trials x xPos x yPos x Time)
%     Binary matrix indicating the spike times for all the trials.
%     The last one is stored in the mosaics of the inner retina object in
%     the responseSpikes slot.
%
% Science and references and computational issues
%
% We set the bipolar scale factor in order to produce a bipolar model that
% generates reasonable RGC spikes for typical viewing conditions.  The sad
% truth is that we don't have a biophysically accurate model of the bipolar
% cells.  Consequently, we have no match for the bipolar current with real
% units.  This is a fudge factor that produces attractive RGC spike rates.
% When we get more information about the bipolar models, we hope to do
% better.
%
% Similarly, the internal calculation converts bipolar current to a
% contrast with a max value of 1.  We can control the max contrast here.
% The code is not set up vary these parameters yet.  We will expose them
% some day.
%
% Example:
%   rgcL.compute(bipolarMosaic??);
%
% See also: rgcMosaic
%
% BW (c) isetbio team

%% Parse inputs
p = inputParser;
p.CaseSensitive = false;

p.addRequired('rgcL',@(x) ~isempty(validatestring(class(x),{'rgcLayer'})));

% For GLM model.  If true, much slower bigger memory
p.addParameter('coupling',false,@islogical);

% Multiple bipolar trials
p.addParameter('bipolarTrials',  [], @(x) isnumeric(x)||iscell(x));
p.addParameter('bipolarScale',50,@isnumeric);
p.addParameter('bipolarContrast',1,@isnumeric);

p.parse(rgcL,varargin{:});
coupling      = p.Results.coupling;
% bipolarTrials = p.Results.bipolarTrials;

% See notes below
bipolarScale    = p.Results.bipolarScale;
bipolarContrast = p.Results.bipolarContrast;

bipolarTrials   = p.Results.bipolarTrials;
if isempty(bipolarTrials)
    bipolarTrials = cell(1,length(rgcL.mosaic));
end
%% Linear stage of the computation

% For now, only deal with one trial case.  Compute the linear response for
% every mosaic. The inputs are already attached.
for ii=1:length(rgcL.mosaic)
    [rgcM,nTrialsLinearResponseM] = ...
        rgcL.mosaic{ii}.computeSeparable(...
        'bipolarContrast',bipolarContrast,...
        'bipolarScale', bipolarScale,...
        'bipolarTrials',bipolarTrials{ii});

    [rgcM,nTrialsSpikeResponseM] = ...
        rgcL.mosaic{ii}.computeSpikes('coupling',coupling,...
        'nTrialsLinear',nTrialsLinearResponseM);
    
    nTrialsSpikes{ii} = nTrialsSpikeResponseM;
end

%% Compute spikes from linear response; possibly for multiple trials

% for ii=1:length(rgcL.mosaic{ii})
%     rgcM = rgcL.mosaic{ii};
%     switch class(rgcM)
%         case {'rgcLinear'}
%             % No linear response implemented yet.
%             disp('No spikes computed for linear RGC mosaic');
%         otherwise
%             % Runs for rgcLNP, rgcGLM
%             % Send the coupling field to decide on the coupling parameter
%             if ~isempty(bipolarTrials)
%                 % Multiple trial case
%                 [rgcL, nTrialsSpikes] = ...
%                     rgcL.computeSpikes(...
%                     'coupling',coupling, ...
%                     'nTrialsLinearResponse',nTrialsLinearResponse);
%             else
%                 % Single trial case
%                 [rgcL, nTrialsSpikes] = rgcL.computeSpikes('coupling',coupling);
%             end
%     end
%
% end
end
