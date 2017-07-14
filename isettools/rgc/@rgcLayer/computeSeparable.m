function [rgcL, nTrialsLinearResponse] = computeSeparable(rgcL, bp, varargin)
% IRCOMPUTELINEARSTSEPARABLE - Computes the RGC mosaic linear response
%
%   ir = irComputeLinearSTSeparable(rgcL, input, 'bipolarTrials', bpTrials)
%
% The linear responses for each mosaic are computed. The linear computation
% is always space-time separation for each mosaic.  The spatial
% computationa, however, is not a convolution because the RF of the cells
% in each mosaic differ.
%
% There are two types of possible input objects.
%
% Required inputs
%   rgcL:    A retina layer object with RGC mosaics
%   bp:      A bipolar cell mosaic object.
%
% Optional inputs
%  bipolarTrials:  Multiple bipolar trials can be sent in using this
%                  variable.
%
% Computational questions
%    * Why do we scale the bipolar voltage input with ieContrast?
%    * Why is the RGC impulse response set to an impulse?  I gather this is
%    because the photocurrent*bipolar equals the observed RGC impulse
%    response?
%
% Computation
%
%  For each mosaic, the center and surround RF responses are calculated
%  (matrix multiply). then the temporal impulse response for the center and
%  surround is calculated.  This continuous operation produces the 'linear'
%  RGC response shown in the window.
%
%  The linear response is the input to irComputeSpikes.
%
% rgcGLM model: The spikes are computed using the recursive influence of
% the post-spike and coupling filters between the nonlinear responses of
% other cells. These computations are carried in irComputeSpikes out using
% code from Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky,
% Simoncelli, Nature, 2008, licensed for modification, which can be found
% at
%
%            http://pillowlab.princeton.edu/code_GLM.html
%
% Outputs: 
%   The linear response and spikes are attached to the inner retina
%   mosaics.  
%
% Examples:
%
%   ir.compute(identityOS);
%
% See also: rgcGLM/rgcCompute, s_vaRGC in WL/WLVernierAcuity
%
% JRG (c) Isetbio team, 2016

%% Check inputs

p = inputParser;
p.CaseSensitive = false;

% ir and bp are both required
p.addRequired('ir',@(x) isequal(class(x),'rgcLayer'));

vFunc = @(x)(isequal(class(x),'bipolarMosaic')||isequal(class(x{1}),'bipolarMosaic'));
p.addRequired('bp',vFunc);

% We can use multiple bipolar trials as input
p.addParameter('bipolarTrials',  [], @(x) isnumeric(x)||iscell(x));
p.addParameter('bipolarScale',  50, @isnumeric);
p.addParameter('bipolarContrast',  1, @isnumeric);

p.parse(rgcL,bp,varargin{:});
bipolarScale = p.Results.bipolarScale;
bipolarContrast = p.Results.bipolarContrast;

%%
bipolarTrials = p.Results.bipolarTrials;
nTrials = 1;
if ~isempty(bipolarTrials) 
    if length(bp)==1
        nTrials = size(bipolarTrials,1); 
    else
        nTrials = size(bipolarTrials{1},1); 
    end
end
        

%% Process the bipolar data
nTrialsLinearResponse = cell(5,1);
for iTrial = 1:nTrials
    % Looping over the rgc mosaics
    for rgcType = 1:length(rgcL.mosaic)                
        % Get the bipolar input data, handle bp mosaic and nTrials cases
        if length(bp) == 1
            if ~isempty(bipolarTrials)
                input   = squeeze(bipolarTrials(iTrial,:,:,:));
            else
                input   = bp.get('response');
            end
        else
            if ~isempty(bipolarTrials)
                input   = squeeze(bipolarTrials{rgcType}(iTrial,:,:,:));
            else
                input   = bp{rgcType}.get('response');
            end
        end
        
        %% Removes the mean and uses a contrast 
        % This is a normalization on the bipolar current.
        % It sc
        % the linear input.  Let's justify or explain or something.
        input = ieContrast(input,'maxC',bipolarContrast);
        
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
        rgcL.mosaic{rgcType}=rgcL.mosaic{rgcType}.set('tCenter all', 1);
        rgcL.mosaic{rgcType}=rgcL.mosaic{rgcType}.set('tSurround all',1);
        
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
    end
end

