function ir = irCompute(ir, outerSegment, varargin)
% Computes the rgc mosaic spikes to an arbitrary stimulus.
%
% The responses for each mosaic are computed one at a time. For a given
% mosaic, first the spatial convolution of the center and surround RFs are
% calculated for each RGB channel, followed by the temporal responses for
% the center and surround and each RGB channel. This results in the linear
% response.
%
% Next, the linear response is put through the generator function. The
% nonlinear response is the input to a function that computes spikes with
% Poisson statistics. For the rgcGLM object, the spikes are computed using
% the recursive influence of the post-spike and coupling filters between
% the nonlinear responses of other cells. These computations are carried
% out using code from Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky,
% Simoncelli, Nature, 2008, licensed for modification, which can be found
% at
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

%%
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) ~isempty(validatestring(class(x),{'ir','irPhys'})));
% validatestring(class(outerSegment),{'osIdentity','osLinear','osBioPhys'})
p.addRequired('outerSegment',@(x) ~isempty(validatestring(class(x),{'osIdentity','osLinear','osBioPhys'})));
p.parse(ir,outerSegment,varargin{:});
osType = class(outerSegment);

% Switch on type of os object
switch osType
    case 'osIdentity'
        %% Identity means straight from the frame buffer to brain   
        % Find properties that haven't been set and set them
        if isempty(osGet(outerSegment,'rgbData'))
            outerSegment = osSet(outerSegment, 'rgbData', rand(64,64,5));
        end        
        if isempty(osGet(outerSegment,'coneSpacing'))
            outerSegment = osSet(outerSegment, 'coneSpacing', 180);
        end        
        if isempty(osGet(outerSegment,'coneSampling'))
            outerSegment = osSet(outerSegment,'coneSampling',.01);
        end
        
        % Run irCreate again?
        
        spTempStim = osGet(outerSegment, 'rgbData');
        
        range = max(spTempStim(:)) - min(spTempStim(:));
        
        % Sometimes the ir class is rgcPhys, which we use for validation.
        % But in general, this is not the case.
        % There is a scientific question about this 10.  We need JRG and EJ
        % to resolve the spiking rate.
        if isequal(class(ir),'irPhys'),   spTempStim = spTempStim./range;
        else                    spTempStim = 10*spTempStim./range;
        end
        
        % Looping over the rgc mosaics
        for rgcType = 1:length(ir.mosaic)
                        
            % We assume a separable space-time receptive field.  This
            % allows us to compute for space first and then time.
            % Space.
            [spResponseCenter, spResponseSurround] = spConvolve(ir.mosaic{rgcType,1}, spTempStim);
            
            % For the subunit model, put each pixel "subunit" of spatial RF
            % through nonlinearity at this point
%             if isa(ir.mosaic{rgcType},'rgcSubunit')
%                 % Change this to get generator function
% %                 spResponseCenter = exp(spResponseCenter);
%                  %                 spResponseSurround = exp(spResponseSurround);
%                  
%                  [sz1,sz2] = size(ir.mosaic{rgcType});
%                  %         obj.spikeResponse{1:sz1,1:sz2,nT,1:nType} = params.value;
%                  
%                  for xc = 1:sz1
%                      for yc = 1:sz2
%                          for nTypeI = 1:nType
%                              obj.responseSpikes{xc,yc,nT+1,nTypeI} = params.value{xc,yc,1,nTypeI};
%                          end
%                      end
%                  end
%             end
            
            % Time. Convolve with the temporal impulse response
            [fullResponse, nlResponse] = ...
                fullConvolve(ir.mosaic{rgcType,1}, spResponseCenter, spResponseSurround);

            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'responseLinear', fullResponse);
            
            % Set the nonlinear response for every rgc subclass except rgcLinear
%             switch class(ir.mosaic{rgcType})
%                 case 'rgcLinear'
%                     % No nonlinear response
%                 otherwise
%                     ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'nlResponse', nlResponse);
%             end
            clear fullResponse nlResponse spResponseCenter spResponseSurround
        end
        
        % Compute spikes
        switch class(ir.mosaic{rgcType})
            case {'rgcLinear','rgcPhys'};
                % No nonlinear response
            otherwise
                for itrial = 1:10
                    ir = irComputeSpikes(ir);
                end
        end       
        
    case {'osLinear'}
        %% Linear OS
        error('Not yet implemented');
    case {'osBioPhys'}
        %% Full biophysical os
        error('Not yet implemented');
        
    otherwise
        error('Unknown os type %s\n',osType);
end


