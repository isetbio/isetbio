function [nTrialsSpikes] = compute(rgcL, varargin)
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
% First,a space-time separable linear response is computed. The computed
% values are stored in the 'responseLinear' slot of each mosaic. The
% critical method is computeSeparable. No noise is added in the linear
% part.
%
% By default, the temporal response is an impulse because of the way the
% bipolar temporal impulse response is set.  Also, the 'dt' of the RGC is
% inherited form the dt of the bipolar calculation.
%
% The spikes are computed from the responseLinear using the computeSpikes
% method. The spiking has a random element.  Running the conversion from
% linear to spikes multiple times produces different spike rasters.
%
% Outputs:
%  nTrialsSpikes:  (trials x xPos x yPos x Time)
%     Binary matrix indicating the spike times for all the trials.
%     The last one is stored in the mosaics of the inner retina object in
%     the responseSpikes slot.
%
% ** Science **
%
% There are a number of unresolved scientific issues about how to represent
% the bipolar responses. We are simply unsure what the right thing to do
% is. Here is what we do.
%
% We convert the bipolar response to a contrast signal.  The maximum value
% of the contrast is set by the parameter 'bipolarContrast', with a default
% of 1.  We compute the inner product between the center and surround
% spatial receptive fields with this contrast.
%
% Then we compute the difference between the center and surround, and we
% scale the output by the parameter 'bipolarScale', which has a default of
% 100.  This value is chosen to influence the spike rates, because this
% setting generates reasonable RGC spikes for typical viewing conditions.
%
% The sad truth is that we don't have a biophysically accurate model of the
% bipolar cells. Consequently, we have no match for the bipolar current
% with real units. This is a HUDGE fudge factor that produces attractive
% RGC spike rates. When we get more information about the bipolar models,
% we hope to do better.
%
% We should be experimenting with bipolarContrast and bipolarScale.  We
% should be reading the literature to try to bring these values into closer
% alignment with the biophysics.  Someone should write the literature.
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
p.addParameter('bipolarTrials',  [], @(x)(isnumeric(x)||iscell(x)));
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
elseif isnumeric(bipolarTrials)
    % Turn the trials into a cell for the compute, below.
    tmp = bipolarTrials;
    clear bipolarTrials
    bipolarTrials{1} = tmp;
    clear tmp
end

%% Linear stage of the computation

% For now, only deal with one trial case.  Compute the linear response for
% every mosaic. The inputs are already attached.
nTrialsSpikes = cell(length(rgcL.mosaic),1);
for ii=1:length(rgcL.mosaic)
    [~,nTrialsLinearResponseM] = ...
        rgcL.mosaic{ii}.computeSeparable(...
        'bipolarContrast',bipolarContrast,...
        'bipolarScale', bipolarScale,...
        'bipolarTrials',bipolarTrials{ii});

    [~,nTrialsSpikeResponseM] = ...
        rgcL.mosaic{ii}.computeSpikes('coupling',coupling,...
        'nTrialsLinear',nTrialsLinearResponseM);
    
    nTrialsSpikes{ii} = nTrialsSpikeResponseM;
end

end
