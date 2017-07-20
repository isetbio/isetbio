function [rgcL, nTrialsLinearResponse] = computeSeparable(rgcM, varargin)
% COMPUTESEPARABLE - Computes the RGC mosaic linear response
%
%   computeSparable(rgcM)
%   @rgcMosaic.computeSeparable(varargin)
%
% Computes the linear responses and spikes for each of the mosaics in a
% retina layer object.
%
% ***
%  PROGRAMMING:  We still have data management to deal with in terms of
%  multiple trials. 
% ***
%
% The linear responses for each mosaic are computed. The linear computation
% is space-time separation for each mosaic.  The spatial computation,
% however, is not a convolution because the RF of the cells within a mosaic
% can differ.
%
% Required inputs
%   rgcM:     A retina mosaic object 
%
% Optional inputs
%  bipolarTrials:  Multiple bipolar trials can be sent in using this
%                  variable.
%
% For each corresponding bp and rgc mosaic, the center and surround RF
% responses are calculated (matrix multiply). then the temporal impulse
% response for the center and surround is calculated.  This continuous
% operation produces the 'linear' RGC response shown in the window.
%
% The spikes are computed from the linear response in a separate routine.
%
% Science and references
%
%    * Why do we scale the bipolar voltage input with ieContrast?
%    * Why is the RGC impulse response set to an impulse?  I gather this is
%    because the photocurrent*bipolar equals the observed RGC impulse
%    response?
%
% rgcGLM model: The spikes are computed using the recursive influence of
% the post-spike and coupling filters between the nonlinear responses of
% other cells. These computations are carried in irComputeSpikes out using
% code from Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky,
% Simoncelli, Nature, 2008, licensed for modification, which can be found
% at
%
%   http://pillowlab.princeton.edu/code_GLM.html
%
% See also: rgcGLM/rgcCompute, s_vaRGC in WL/WLVernierAcuity
%
% JRG/BW (c) Isetbio team, 2016

%% Check inputs

p = inputParser;
p.CaseSensitive = false;

% ir and bp are both required
p.addRequired('rgcM',@(x) isequal(class(x),'rgcMosaic'));

% We can use multiple bipolar trials as input
p.addParameter('bipolarScale',  50, @isnumeric);
p.addParameter('bipolarContrast',  1, @isnumeric);

% Not used now.  To be added later.
p.addParameter('bipolarTrials',  [], @(x) isnumeric(x)||iscell(x));

p.parse(rgcM,varargin{:});
bipolarScale    = p.Results.bipolarScale;
bipolarContrast = p.Results.bipolarContrast;

%%
% bipolarTrials = p.Results.bipolarTrials;
% nTrials = 1;
% if ~isempty(bipolarTrials)
%     if length(bpMosaic)==1
%         nTrials = size(bipolarTrials,1);
%     else
%         nTrials = size(bipolarTrials{1},1);
%     end
% end


%% Process the bipolar data
% nTrialsLinearResponse = cell(5,1);
% for iTrial = 1:nTrials
% Looping over the rgc mosaics
% for rgcType = 1:length(rgcL.mosaic)
%     % Get the bipolar input data, handle bp mosaic and nTrials cases
%     if length(bpMosaic) == 1
%         if ~isempty(bipolarTrials)
%             input   = squeeze(bipolarTrials(iTrial,:,:,:));
%         else
%             input   = bpMosaic.get('response');
%         end
%     else
%         if ~isempty(bipolarTrials)
%             input   = squeeze(bipolarTrials{rgcType}(iTrial,:,:,:));
%         else
%             input   = bpMosaic{rgcType}.get('response');
%         end
%     end

%% Removes the mean and uses a contrast
% This is a normalization on the bipolar current.
% It sc
% the linear input.  Let's justify or explain or something.

input = ieContrast(rgcM.get('response'),'maxC',bipolarContrast);

% vcNewGraphWin; ieMovie(stim);
% foo = zeros(10,11,size(stim,3));
% for ii=1:size(stim,3), foo(:,:,ii) = imresize(stim(:,:,ii),[10,11]); end
% ieMovie(foo);

%% Set the rgc impulse response to an impulse
% When we feed a bipolar object into the inner retina, we don't
% need to do temporal convolution. We have the tCenter and
% tSurround properties for the rgcMosaic, so we set them to an
% impulse to remind us that the temporal repsonse is already
% computed.
rgcL.mosaic{rgcType} = rgcL.mosaic{rgcType}.set('tCenter all', 1);
rgcL.mosaic{rgcType} = rgcL.mosaic{rgcType}.set('tSurround all',1);

% We use a separable space-time receptive field that computes for
% space here.  We will implement temporal response later.
[respC, respS] = rgcSpaceDot(rgcL.mosaic{rgcType}, input);
% vcNewGraphWin; ieMovie(respC);

%% Convolve with the temporal impulse response
% If the temporal IR is an impulse, we don't need to do the
% temporal convolution.  I'm leaving this here as a reminder that
% we need to do this if we run any of EJ's original GLM models
% (image -> spikes with no cone mosaic or bipolar).
% if ir.mosaic{rgcType}.tCenter{1,1} ~= 1
%     respC = timeConvolve(ir.mosaic{rgcType}, respC, 'c');
%     respS = timeConvolve(ir.mosaic{rgcType}, respS, 's');
% end
% ieMovie(respC - respS);


%% Deal with multiple trial issues
if ~isempty(bipolarTrials)
    if iTrial == 1
        nTrialsLinearResponse{rgcType} = zeros([nTrials,size(respC)]);
    end
    nTrialsLinearResponse{rgcType}(iTrial,:,:,:) =  bipolarScale*(respC - respS);
end

% Store the last trial
if iTrial == nTrials
    % Store the linear response
    rgcL.mosaic{rgcType} = mosaicSet(rgcL.mosaic{rgcType},'response linear', bipolarScale*(respC - respS));
    % rgcL.mosaic{rgcType}.window;
end
% end
end


