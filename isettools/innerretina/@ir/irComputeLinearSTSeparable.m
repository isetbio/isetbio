function [ir, nTrialsLinearResponse] = irComputeLinearSTSeparable(ir, bp, varargin)
% Computes the mosaic's linear response to an input
%
%   ir = irComputeLinearSTSeparable(ir, input, varargin)
%
% The linear responses for each mosaic are computed one at a time. The
% linear computation is always space-time separable in here?
%
% There are two types of possible input objects.
%
%      'osDisplayRGB' - frame buffer values for a spatiotemporal stimulus
%           stored in an outer segment object.
%      'bipolar' - the bipolar cell object with a signal that has gone
%           through temporal filtering and possibly spatial subunit
%           nonlinearities.
% 
% For a given mosaic, first the spatial convolution of the center and
% surround RFs are calculated for each RGB channel, followed by the
% temporal responses for the center and surround and each RGB channel. This
% results in the linear response.
%
% The linear response is the input to irComputeSpikes.
%
% For the rgcGLM object, the spikes can be computed using the recursive
% influence of the post-spike and coupling filters between the nonlinear
% responses of other cells. These computations are carried in
% irComputeSpikes out using code from Pillow, Shlens, Paninski, Sher,
% Litke, Chichilnisky, Simoncelli, Nature, 2008, licensed for modification,
% which can be found at
%
% http://pillowlab.princeton.edu/code_GLM.html
%
% Outline:
%  * Normalize stimulus
%  * Compute linear response
%     - spatial convolution
%     - temporal convolution
%  * Compute nonlinear response
% [spiking responses are calculated by subclass versions of rgcCompute]
%
% Inputs: inner retina object, outersegment object.
%
% Outputs: the inner object with fixed linear and noisy spiking responses.
%
% Example:
%   ir.compute(identityOS);
%   irCompute(ir,identityOS);
%
% See also: rgcGLM/rgcCompute
%
% JRG (c) isetbio team

%% Parse
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) isequal(class(x),'ir')||isequal(class(x),'irPhys'));
vFunc = @(x)(isequal(class(x),'bipolar')||isequal(class(x{1}),'bipolar'));
p.addRequired('bp',vFunc);

p.addParameter('bipolarTrials',  [], @isnumeric);

p.parse(ir,bp,varargin{:});

ir = p.Results.ir;
bp = p.Results.bp;

bipolarTrials = p.Results.bipolarTrials;

%% Get the input data

% Possible osTypes are osIdentity, osLinear, and osBiophys
% Only osIdentity is implemented now.
if length(bp) == 1
    osType = class(bp);
else
    osType = class(bp{1});
end
% Bipolar test case in t_coneMosaic
% t_rgcBar, others to be named.

if ~isempty(bipolarTrials)
    nTrials = size(bipolarTrials,1);
else
    nTrials = 1;
end

for iTrial = 1:nTrials
    
    % Looping over the rgc mosaics
    for rgcType = 1:length(ir.mosaic)
        
        % Determine the range of the rgb input data
        if length(bp) == 1
            stim   = bp.get('response');
        else
            stim   = bp{rgcType}.get('responseCenter');
        end
        switch class(ir.mosaic{rgcType})
            case 'rgcPhys'
                magFactor = 7.9; % due to bipolar filter
                stim = magFactor*ieContrast(stim);
            otherwise
                stim = ieContrast(stim);
        end
        % ieMovie(stim);
        
        % Set the rgc impulse response to an impulse
        ir.mosaic{rgcType}=ir.mosaic{rgcType}.set('tCenter all', 1);
        ir.mosaic{rgcType}=ir.mosaic{rgcType}.set('tSurround all',0);
        
        % We use a separable space-time receptive field.  This allows
        % us to compute for space first and then time. Space.
        [respC, respS] = spConvolve(ir.mosaic{rgcType}, stim);
        % ieMovie(respC);
        
        % Convolve with the temporal impulse response
        respC = timeConvolve(ir.mosaic{rgcType}, respC, 'c');
        respS = timeConvolve(ir.mosaic{rgcType}, respS, 's');
        % Delete fullConvolve
        % ieMovie(respC - respS);
        
        
        if ~isempty(bipolarTrials)
            if iTrial == 1
                nTrialsLinearResponse = zeros([nTrials,size(respC)]);
            end
            nTrialsLinearResponse(iTrial,:,:,:) =  respC - respS;
        end
        
        if iTrial == nTrials
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'response linear', respC - respS);
        end
        
    end
end



