function [ir, nTrialsLinearResponse] = irComputeLinearSTSeparable(ir, input, varargin)
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

% allowableInputs = {'osDisplayRGB','bipolar'};
% p.addRequired('input',@(x) ismember(class(x),allowableInputs));
if length(input) == 1
    vFunc = @(x) ~isempty(validatestring(class(x),{'osDisplayRGB','bipolar'}));
else
    vFunc = @(x) ~isempty(validatestring(class(x{1}),{'osDisplayRGB','bipolar'}));
end
p.addRequired('input',vFunc);

p.addParameter('bipolarTrials',  [], @isnumeric);

p.parse(ir,input,varargin{:});

%
ir    = p.Results.ir;
input = p.Results.input;

bipolarTrials = p.Results.bipolarTrials;

%% Get the input data

% Possible osTypes are osIdentity, osLinear, and osBiophys
% Only osIdentity is implemented now.
if length(input) == 1
    osType = class(input);
else
    osType = class(input{1});
end
switch osType
    case 'osDisplayRGB'
        % Display RGB means straight from the frame buffer to your brain
        % @JRG:  Test cases are in EJ Figure reproduction
        % Others to be named!
        
        spTempStim = osGet(input, 'rgb data');
        if isempty(spTempStim), error('No rgb data'); end
        
        % The Pillow code expects the input to be normalized to the mean of
        % zero, like a contrast. In this way a frame buffer of 1024 or 512
        % or 256 would all get scaled into the same [-0.5, 0.5] range.
        stim = ieContrast(spTempStim);
        %         range = max(spTempStim(:)) - min(spTempStim(:));
        %         spTempStim = (spTempStim - mean(spTempStim(:)))/range;
        
        % Looping over each rgc
        for rgcType = 1:length(ir.mosaic)                  
            
            % We use a separable space-time receptive field.  This allows
            % us to compute for space first and then time. Space.
            [respC, respS] = spConvolve(ir.mosaic{rgcType}, stim);
            % ieMovie(respC);
            
            tResponse = ir.mosaic{rgcType}.tCenter{1,1};
            
            ir.mosaic{rgcType}=ir.mosaic{rgcType}.set('tCenter all', tResponse);
            ir.mosaic{rgcType}=ir.mosaic{rgcType}.set('tSurround all',0);        
            
            % Convolve with the temporal impulse response
            respC = timeConvolve(ir.mosaic{rgcType}, respC, 'c');
            respS = timeConvolve(ir.mosaic{rgcType}, respS, 's');
            % Delete fullConvolve
            % ieMovie(respC - respS);
                        
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'response linear', respC - respS);
            
        end       
        
    case {'bipolar'}
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
                if length(input) == 1
                    stim   = bipolarGet(input, 'response');
                else
                    stim   = bipolarGet(input{rgcType}, 'response');
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
    otherwise
        error('Unknown os type %s\n',osType);
end


