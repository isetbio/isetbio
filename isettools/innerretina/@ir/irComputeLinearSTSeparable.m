function [ir, nTrialsLinearResponse] = irComputeLinearSTSeparable(ir, bp, varargin)
% IRCOMPUTELINEARSTSEPARABLE - Computes the RGC mosaic linear response
%
%   ir = irComputeLinearSTSeparable(ir, input, 'bipolarTrials', bpTrials)
%
% The linear responses for each mosaic are computed. The linear computation
% is always space-time separation for each mosaic.  The spatial
% computationa, however, is not a convolution because the RF of the cells
% in each mosaic differ.
%
% There are two types of possible input objects.
%
% Required inputs
%   ir:      An inner retina object with RGC mosaics
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
%   rgcGLM model: The spikes are computed using the recursive influence of
%   the post-spike and coupling filters between the nonlinear responses of
%   other cells. These computations are carried in irComputeSpikes out
%   using code from Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky,
%   Simoncelli, Nature, 2008, licensed for modification,
%   which can be found at
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
p.addRequired('ir',@(x) isequal(class(x),'ir')||isequal(class(x),'irPhys'));

vFunc = @(x)(isequal(class(x),'bipolar')||isequal(class(x{1}),'bipolar'));
p.addRequired('bp',vFunc);

% We can use multiple bipolar trials as input
p.addParameter('bipolarTrials',  [], @isnumeric);

p.parse(ir,bp,varargin{:});

bipolarTrials = p.Results.bipolarTrials;
nTrials = 1;
if ~isempty(bipolarTrials), nTrials = size(bipolarTrials,1); end

%% Process the bipolar data

for iTrial = 1:nTrials
    
    % Looping over the rgc mosaics
    for rgcType = 1:length(ir.mosaic)
        
        % BW inserted.  I think these are supposed to match.
        ir.mosaic{rgcType}.dt = bp.timeStep;
        
        % Determine the range of the rgb input data
        if length(bp) == 1,   stim   = bp.get('response');
        else,                 stim   = bp{rgcType}.get('responseCenter');
        end
        
        % JRG removes the mean and uses a contrast (or scaled contrast) as
        % the linear input.  Let's justify or explain or something.
        switch class(ir.mosaic{rgcType})
            case 'rgcPhys'
                % Let's get rid of this parameter.  Or explain.  Or
                % something.
                magFactor = 7.9; % due to bipolar filter
                stim = magFactor*ieContrast(stim);
            otherwise
                stim = ieContrast(stim);
        end
        % ieMovie(stim);
        
        % Set the rgc impulse response to an impulse
        % Why? (BW)  And if this is right, then why would we even apply it?
        ir.mosaic{rgcType}=ir.mosaic{rgcType}.set('tCenter all', 1);
        ir.mosaic{rgcType}=ir.mosaic{rgcType}.set('tSurround all',0);
        
        % We use a separable space-time receptive field.  This allows
        % us to compute for space first and then time. Space.  But
        % spConvolve is not a convolution.  So, we should rename this
        % routine to just rgcSpace(), I think (BW).
        [respC, respS] = spConvolve(ir.mosaic{rgcType}, stim);
        % ieMovie(respC);
        
        % Convolve with the temporal impulse response
        % Why, given that these are impulses?  What's going on here?
        respC = timeConvolve(ir.mosaic{rgcType}, respC, 'c');
        respS = timeConvolve(ir.mosaic{rgcType}, respS, 's');
        % ieMovie(respC - respS);
        
        % Deal with multiple trial issues
        if ~isempty(bipolarTrials)
            if iTrial == 1
                nTrialsLinearResponse = zeros([nTrials,size(respC)]);
            end
            nTrialsLinearResponse(iTrial,:,:,:) =  respC - respS;
        end
        
        % Store the last trial
        if iTrial == nTrials
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'response linear', respC - respS);
        end
    end
end

