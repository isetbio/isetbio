function [rgcM, nTrialsLinearResponse] = computeSeparable(rgcM, varargin)
% COMPUTESEPARABLE - Computes the RGC mosaic linear response
%
%   computeSparable(rgcM)
%   @rgcMosaic.computeSeparable(varargin)
%
% Computes the linear responses and spikes for each of the mosaics in a
% retina layer object.
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
%                  variable. It is a cell array {nTrials,x,y,t}
%
% This methods produces the 'linear' RGC response shown in the window.
%
% For each corresponding bipolarMosaic and rgcMosaic, the center and
% surround RF responses are calculated (matrix multiply). In principle, the
% temporal impulse response for the center and surround is calculated but
% in practice for now we only use the impulse.  This maintains the
% bipolar/cone impulse response to match the RGC impulse response. 
%
% The spikes are computed from these linear responses in the computeSpikes
% method.
%
% Science and references - NEEDS TO BE EDITED (BW)
%
% See @rgcLayer.compute for a discussion
%
%    * Why do we scale the bipolar voltage input with ieContrast?
%    * Why is the RGC impulse response set to an impulse?  Because the
%    photocurrent*bipolar equals the observed RGC impulse response?
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

validTypes = {'rgcGLM','rgcLNP'};
vFunc = @(x)(ismember(class(x),validTypes));
p.addRequired('rgcM',vFunc);

% We can use multiple bipolar trials as input
p.addParameter('bipolarScale',  50, @isnumeric);
p.addParameter('bipolarContrast',  1, @isnumeric);
p.addParameter('bipolarTrials',  [], @(x) isnumeric(x)||iscell(x));

p.parse(rgcM,varargin{:});
bipolarScale    = p.Results.bipolarScale;
bipolarContrast = p.Results.bipolarContrast;
bipolarTrials = p.Results.bipolarTrials;

%% Loop over multiple trials

nTrials = 1;
if ~isempty(bipolarTrials), nTrials = size(bipolarTrials,1); end

for iTrial = 1:nTrials
    %% Removes the mean of the bipolar mosaic input, converts to contrast
    
    % This is a normalization on the bipolar current.
    % Let's justify or explain or something.
    
    if ~isempty(bipolarTrials)
        % Set the contrast for this trial
        input   = ieContrast(squeeze(bipolarTrials(iTrial,:,:,:)),'maxC',bipolarContrast);
    else
        input   = ieContrast(rgcM.input.get('response'),'maxC',bipolarContrast);
    end
    % vcNewGraphWin; ieMovie(input);
    
    %% Set the rgc impulse response to an impulse
    
    % When we feed a bipolar object into the inner retina, we don't need to
    % do temporal convolution. We have the tCenter and tSurround properties
    % for the rgcMosaic, so we set them to an impulse to remind us that the
    % temporal repsonse is already computed.
    rgcM.set('tCenter all', 1);
    rgcM.set('tSurround all',1);
    
    % We use a separable space-time receptive field that computes for
    % space here.  We will implement temporal response later.
    [respC, respS] = rgcSpaceDot(rgcM, input);
    % vcNewGraphWin; ieMovie(respC);
    
    % Convolve with the temporal impulse response
    %
    % We commented out for now, but this temporal convolution should be
    % allowed at some point.
    %
    % If the temporal IR is an impulse, we don't need to do the
    % temporal convolution.  I'm leaving this here as a reminder that
    % we need to do this if we run any of EJ's original GLM models
    % (image -> spikes with no cone mosaic or bipolar).
    % if ir.mosaic{rgcType}.tCenter{1,1} ~= 1
    %     respC = timeConvolve(ir.mosaic{rgcType}, respC, 'c');
    %     respS = timeConvolve(ir.mosaic{rgcType}, respS, 's');
    % end
    % ieMovie(respC - respS);
    
    
    %% Deal with multiple trial returns
    
    if isempty(bipolarTrials)
        % Shouldn't this be if nTrial == 1?
        % Then if iTrial == 1 and nTrial > 1 we allocate?
        nTrialsLinearResponse = [];
        
    elseif ~isempty(bipolarTrials) && iTrial == 1
        [nr,nc,nt] = size(respC);
        % Allocate
        nTrialsLinearResponse  =  zeros(nTrials, nr, nc, nt);
        % Put in the first trial
        nTrialsLinearResponse(iTrial,:,:,:) =  bipolarScale*(respC - respS);
        
    else
        % I don't understand why Matlab thinks this is incrementing on the
        % loop. (BW)
        nTrialsLinearResponse(iTrial,:,:,:) =  bipolarScale*(respC - respS);
    end
    
    % Store the last trial
    if iTrial == nTrials
        % Store the last linear response in the object
        rgcM.set('response linear', bipolarScale*(respC - respS));
    end
end

end


