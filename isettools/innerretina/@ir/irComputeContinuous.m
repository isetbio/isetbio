function ir = irComputeContinuous(ir, outerSegment, varargin)
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

p.addRequired('ir',@(x) isequal(class(x),'ir'));
% validatestring(class(outerSegment),{'osIdentity','osLinear','osBioPhys'})
p.addRequired('outerSegment',@(x) ~isempty(validatestring(class(x),{'osIdentity','osLinear','osBioPhys'})));
% p.parse(ir,outerSegment,varargin{:});
osType = class(outerSegment);
%% Get the input data

% Possible osTypes are osIdentity, osLinear, and osBiophys
% Only osIdentity is implemented now.
osType = class(outerSegment);
switch osType
    case 'osDisplayRGB'
        %% Identity means straight from the frame buffer to brain
        % Find properties that haven't been set and set them
        if isempty(osGet(outerSegment,'rgbData'))
            outerSegment = osSet(outerSegment, 'rgbData', rand(64,64,5));
        end
        if isempty(osGet(outerSegment,'patchSize'))
            outerSegment = osSet(outerSegment, 'patchSize', 180);
        end
        if isempty(osGet(outerSegment,'timeStep'))
            outerSegment = osSet(outerSegment,'timeStep',.01);
        end
        
        
        %% Linear computation       
        
        % Determine the range of the rgb input data
        spTempStim = osGet(outerSegment, 'rgbData');
        range = max(spTempStim(:)) - min(spTempStim(:));
        
        % Special case. If the ir class is rgcPhys, which we use for
        % validation. But in general, this is not the case. There is a
        % scientific question about this 10.  We need JRG and EJ to resolve
        % the spiking rate.
        % James needs to change the spatial RF and temporal weights in order to
        % make these models have the right spike rate, and the 10 is a hack to
        % approximate that.
        if isequal(class(ir),'irPhys'),   spTempStim = spTempStim./range - mean(spTempStim(:))/range;
        else                    spTempStim = 300*(spTempStim./range - mean(spTempStim(:))/range);
%             else                    spTempStim = 10*(spTempStim./range);
        end
        
        % Looping over the rgc mosaics
        for rgcType = 1:length(ir.mosaic)
            
            % We use a separable space-time receptive field.  This allows
            % us to compute for space first and then time. Space.
            [spResponseCenter, spResponseSurround] = spConvolve(ir.mosaic{rgcType,1}, spTempStim);
            
            % For the subunit model, put each pixel "subunit" of spatial RF
            % through nonlinearity at this point
%             if isa(ir.mosaic{rgcType},'rgcSubunit')
%                 % Change this to get generator function
%                 
%                 modeltype = 'pixel';
%                 % modeltype = 'surround';
%                 [spResponseCenter, spResponseSurround] = ...
%                     subunitPooling(spResponseCenter, spResponseSurround, 'model', modeltype);
%                 
%             end
            
            % Convolve with the temporal impulse response
            responseLinear = ...
                fullConvolve(ir.mosaic{rgcType,1}, spResponseCenter, spResponseSurround);
            
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'responseLinear', responseLinear);
            
        end
        
    case {'osLinear','osBioPhys'}
        %% Linear OS
        
        % Programming TODO: Connect L, M and S cones to RGC centers and
        % surrounds appropriately
               
        % Determine the range of the cone current
        spTempStim = osGet(outerSegment, 'coneCurrentSignal');
        
        % Probably need to do this by cone type
        range = max(spTempStim(:)) - min(spTempStim(:));
        
        % Special case. If the ir class is rgcPhys, which we use for
        % validation. But in general, this is not the case. There is a
        % scientific question about this 10.  We need JRG and EJ to resolve
        % the spiking rate.
        % James needs to change the spatial RF and temporal weights in order to
        % make these models have the right spike rate, and the 10 is a hack to
        % approximate that.
        spTempStim = spTempStim./range;
        spTempStim = spTempStim - mean(spTempStim(:));
        
        % Looping over the rgc mosaics
        for rgcType = 1:length(ir.mosaic)
            
            % We use a separable space-time receptive field.  This allows
            % us to compute for space first and then time. Space.
            [spResponseCenter, spResponseSurround] = spConvolve(ir.mosaic{rgcType,1}, spTempStim);
            
            % For the subunit model, put each pixel "subunit" of spatial RF
            % through nonlinearity at this point
            if isa(ir.mosaic{rgcType},'rgcSubunit')
                % Change this to get generator function
                spResponseCenter = cellfun(@exp,spResponseCenter,'uniformoutput',false);
                spResponseSurround = cellfun(@exp,spResponseSurround,'uniformoutput',false);
            end
            
            % Convolve with the temporal impulse response
            % responseLinear = ...
            %    fullConvolve(ir.mosaic{rgcType,1}, spResponseCenter, spResponseSurround);
            
            % spResponseSum = cellfun(@plus,spResponseCenter ,spResponseSurround,'uniformoutput',false);
            % spResponseVecRS = cellfun(@sum,(cellfun(@sum,spResponseSum,'uniformoutput',false)),'uniformoutput',false);
            % spResponseVec=cellfun(@squeeze,spResponseVecRS,'uniformoutput',false);
            
            cellCtr = 0;
            nCells = size(ir.mosaic{rgcType}.cellLocation);
            for xc = 1:nCells(1)
                for yc = 1:nCells(2)
                    spResponseFull = spResponseCenter{xc,yc} + spResponseSurround{xc,yc};
                    spResponseVec{xc,yc} = squeeze(mean(mean(spResponseFull,1),2))';
                end
            end
            
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'responseLinear', spResponseVec);
            
        end
%     case {'osBioPhys'}
%         %% Full biophysical os
%         error('Not yet implemented');

    case {'bipolar'}

        %% Linear computation       
        
        % Determine the range of the rgb input data
        spTempStimCenter = bipolarGet(outerSegment, 'responseCenter');        
        spTempStimSurround = bipolarGet(outerSegment, 'responseSurround');
        
        rangeCenter = max(spTempStimCenter(:)) - min(spTempStimCenter(:));
        rangeSurround = max(spTempStimSurround(:)) - min(spTempStimSurround(:));
        
        spTempStimCenter = spTempStimCenter./rangeCenter - mean(spTempStimCenter(:))/rangeCenter;
        spTempStimSurround = spTempStimSurround./rangeSurround - mean(spTempStimSurround(:))/rangeSurround;
        
        % Looping over the rgc mosaics
        for rgcType = 1:length(ir.mosaic)
            
            % We use a separable space-time receptive field.  This allows
            % us to compute for space first and then time. Space.
            [spResponseCenter, spResponseSurround] = spConvolve(ir.mosaic{rgcType,1}, spTempStimCenter, spTempStimSurround);           
            
            
            szCenter = size(spResponseCenter);
            for s1 = 1:szCenter(1)
                for s2 = 1:szCenter(2)
%                     spResponseCenter{s1,s2}(isnan(spResponseCenter{s1,s2})) = 0;
%                     spResponseSurround{s1,s2}(isnan(spResponseSurround{s1,s2})) = 0;
                    spResponseCenter{s1,s2} = 40*spResponseCenter{s1,s2};
                    spResponseSurround{s1,s2} = 40*spResponseSurround{s1,s2};
                end
            end
%             
%             figure;
%             hold on;
%             for s1 = 1:szCenter(1)
%                 for s2 =1:szCenter(2)
% %                     plot(squeeze(spResponseCenter{s1,s2}(1,1,:,1)));
% %                     plot(squeeze(spResponseSurround{s1,s2}(1,1,:,1)));
%                     mx(s1,s2) = max(squeeze(spResponseCenter{s1,s2}(1,1,:,1))+squeeze(spResponseSurround{s1,s2}(1,1,:,1)));
% %                     plot((squeeze(spResponseCenter{s1,s2}(1,1,:,1))+squeeze(spResponseSurround{s1,s2}(1,1,:,1)))./mx(s1,s2));
% %                     spResponseCenter{s1,s2} = 10000*spResponseCenter{s1,s2}./mx(s1,s2);
% %                     spResponseSurround{s1,s2} = 10000*spResponseSurround{s1,s2}./mx(s1,s2);
%                     plot((squeeze(spResponseCenter{s1,s2}(1,1,:,1))) + squeeze(spResponseSurround{s1,s2}(1,1,:,1)));
% %                     plot((squeeze(spResponseCenter{s1,s2}(1,1,:,1))));
% %                     hold on; plot(squeeze(spResponseSurround{s1,s2}(1,1,:,1)));
%                 end
%             end
            
            % Convolve with the temporal impulse response
            responseLinear = ...
                fullConvolve(ir.mosaic{rgcType,1}, spResponseCenter, spResponseSurround);
            
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'responseLinear', responseLinear);


%             figure; 
%             hold on;
%             for s1 = 1:szCenter(1)
%                 for s2 = 1:szCenter(2)
% %                     plot(squeeze(spResponseCenter{s1,s2}(1,1,:,1)));
% %                     plot(squeeze(spResponseSurround{s1,s2}(1,1,:,1)));
%                         plot(squeeze(responseLinear{s1,s2,2}(1,1,:)));
%                 end
%             end
            
%             xlabel('time (msec)','fontsize',14); ylabel('Activation','fontsize',14);
%             title('Linear Activation before temporal filtering');
        end
    otherwise
        error('Unknown os type %s\n',osType);
end





