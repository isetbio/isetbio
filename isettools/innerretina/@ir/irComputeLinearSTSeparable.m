function ir = irComputeLinearSTSeparable(ir, input, varargin)
% Computes the mosaic's linear response to an input
%
%   ir = irComputeLinearSTSeparable(ir, input, varargin)
%
% The linear responses for each mosaic are computed one at a time. The
% linear computation is always space-time separable in here?
%
% See irCompute for a listing of the possible input types.
% 
% For a
% given mosaic, first the spatial convolution of the center and surround
% RFs are calculated for each RGB channel, followed by the temporal
% responses for the center and surround and each RGB channel. This results
% in the linear response.
%
% Next, the linear response is put through the generator function. The
% nonlinear response is the input to a function that computes spikes with
% Poisson statistics. 
%
% For the rgcGLM object, the spikes are computed using the recursive
% influence of the post-spike and coupling filters between the nonlinear
% responses of other cells. These computations are carried out using code
% from Pillow, Shlens, Paninski, Sher, Litke, Chichilnisky, Simoncelli,
% Nature, 2008, licensed for modification, which can be found at
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

p.addRequired('ir',@(x) isequal(class(x),'ir'));

allowableInputs = {'osDisplayRGB','bipolar'};
p.addRequired('input',@(x) ismember(class(x),allowableInputs));

p.parse(ir,input,varargin{:});
ir = p.Results.ir;
input = p.Results.input;

%% Get the input data

% Possible osTypes are osIdentity, osLinear, and osBiophys
% Only osIdentity is implemented now.
osType = class(input);
switch osType
    case 'osDisplayRGB'
        % Display RGB means straight from the frame buffer to your brain
        % @JRG:  Test cases are in EJ Figure reproduction
        % Others to be named!
        
        spTempStim = osGet(input, 'rgbData');
        if isempty(spTempStim), error('No rgb data'); end
        
        % The Pillow code expects the input to be normalized to the mean of
        % zero, like a contrast. In this way a frame buffer of 1024 or 512
        % or 256 would all get scaled into the same [-0.5, 0.5] range.
        spTempStim = ieContrast(spTempStim);
        %         range = max(spTempStim(:)) - min(spTempStim(:));
        %         spTempStim = (spTempStim - mean(spTempStim(:)))/range;
        
        % Looping over each rgc
        for rgcType = 1:length(ir.mosaic)
            
            % We use a separable space-time receptive field.  This allows
            % us to compute for space first and then time. Space.
            [spResponseCenter, spResponseSurround] = ...
                spConvolve(ir.mosaic{rgcType,1}, spTempStim);
            
            % Convolve with the temporal impulse response
            responseLinear = ...
                fullConvolve(ir.mosaic{rgcType,1}, spResponseCenter, spResponseSurround);
            
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'response linear', responseLinear);
            
        end
        
        %     case {'osLinear','osBioPhys'}
        %         % Linear OS
        %         % Test scripts for this case could be
        %         %
        %
        %         % Programming TODO: Connect L, M and S cones to RGC centers and
        %         % surrounds appropriately
        %
        %         % Determine the range of the cone current
        %         spTempStim = osGet(input, 'cone current signal');
        %
        %         % Probably need to do this by cone type
        %         % We should make this a function with a meaningful name
        %         % We have ieScale.  But maybe there should be ieContrast
        %         range = max(spTempStim(:)) - min(spTempStim(:));
        %         spTempStim = (spTempStim - mean(spTempStim(:)))/range;
        %
        %         % Looping over the rgc mosaics
        %         for rgcType = 1:length(ir.mosaic)
        %
        %             % We use a separable space-time receptive field.  This allows
        %             % us to compute for space first and then time. Space.
        %             [spResponseCenter, spResponseSurround] = spConvolve(ir.mosaic{rgcType,1}, spTempStim);
        %
        %             spResponseVec = cell(nCells(1),nCells(2));
        %             nCells = size(ir.mosaic{rgcType}.cellLocation);
        %             for xc = 1:nCells(1)
        %                 for yc = 1:nCells(2)
        %                     spResponseFull = spResponseCenter{xc,yc} + spResponseSurround{xc,yc};
        %                     spResponseVec{xc,yc} = squeeze(mean(mean(spResponseFull,1),2))';
        %                 end
        %             end
        %
        %             % Store the linear response
        %             ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'responseLinear', spResponseVec);
        %
        %         end
        %         %     case {'osBioPhys'}
        %         %         %% Full biophysical os
        %         %         error('Not yet implemented');
        
    case {'bipolar'}
        % Bipolar test case in t_coneMosaic
        % t_rgcBar, others to be named.
        
        % Determine the range of the rgb input data
        spTempStimCenter   = bipolarGet(input, 'response center');
        spTempStimCenter   = ieContrast(spTempStimCenter);
        spTempStimSurround = bipolarGet(input, 'response surround');
        spTempStimSurround = ieContrast(spTempStimSurround);
        
        %         % This will be the ieContrast command
        %         rangeCenter = max(spTempStimCenter(:))     - min(spTempStimCenter(:));
        %         rangeSurround = max(spTempStimSurround(:)) - min(spTempStimSurround(:));
        %
        %         if rangeCenter ~= 0
        %             spTempStimCenter = spTempStimCenter./rangeCenter - mean(spTempStimCenter(:))/rangeCenter;
        %         end
        %         if rangeSurround ~= 0
        %             spTempStimSurround = spTempStimSurround./rangeSurround - mean(spTempStimSurround(:))/rangeSurround;
        %         end
        %
        
        % Looping over the rgc mosaics
        for rgcType = 1:length(ir.mosaic)
                
            % Set the rgc impulse response to an impulse
            % The reason for a negative one instead of a positive one is
            % @JRG to explain
            ir.mosaic{rgcType}.set('tCenter all', -1);
            ir.mosaic{rgcType}.set('tSurround all',0);
            
            %             xNum = size(ir.mosaic{rgcType}.sRFcenter,1);
            %             yNum = size(ir.mosaic{rgcType}.sRFcenter,2);
            %
            %             tCenterBP    = -1;
            %             tCenterNew   = cell(xNum,yNum);
            %             tSurroundNew = cell(xNum,yNum);
            %             for xcell = 1:xNum
            %                 for ycell = 1:yNum
            %                     tCenterNew{xcell,ycell} = tCenterBP;
            %                     tSurroundNew{xcell,ycell} = 0;
            %                 end
            %             end
            %
            %             ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'tCenter',tCenterNew);
            %             ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'tSurround',tSurroundNew);
            
            % @JRG may need for the EJ repository
            %
            % szspt = size(spTempStimCenter);
            % MAY NEED TO ALLOW RESIZE FOR t_rgcCascade, t_rgcPeriphery
            %                 if 0 %isa(ir.mosaic{rgcType},'rgcPhys') && (szspt(1) ~= 80 && szspt(1) ~= 40)
            %                      %                 spTempStimCenterRS = spTempStimCenter;
            %                      %                 spTempStimSurroundRS = spTempStimSurround;
            %                     spTempStimCenterRS = zeros(80,40,size(spTempStimCenter,3));
            %                     spTempStimSurroundRS = zeros(80,40,size(spTempStimCenter,3));
            %                     for frstim = 1:size(spTempStimCenter,3)
            %                         spTempStimCenterRS(:,:,frstim) =    imresize(spTempStimCenter  (:,:,frstim),[80 40],'method','nearest');
            %                         spTempStimSurroundRS(:,:,frstim) =  imresize(spTempStimSurround(:,:,frstim),[80 40],'method','nearest');
            %                     end
            %
            %                     clear spTempStimCenter spTempStimSurround
            %
            %                     spTempStimCenter = spTempStimCenterRS;
            %                     spTempStimSurround = spTempStimSurroundRS;
            %                 end
            
            
            
            % We use a separable space-time receptive field.  This allows
            % us to compute for space first and then time. Space.
            [spResponseCenter, spResponseSurround] = ...
                spConvolve(ir.mosaic{rgcType,1}, spTempStimCenter, spTempStimSurround);
            
            
            %             for cellNum = 1:length(ir.mosaic{rgcType}.sRFcenter)
            %                 linEqDiscBig(:,:,frstim) =    imresize(linEqDisc{cellNum,1}(:,:,frstim),[80 40]);
            %             end
            
            % spResponseSurround{1,1} = zeros(size(spResponseCenter{1,1}));
            % Convolve with the temporal impulse response
            responseLinear = ...
                fullConvolve(ir.mosaic{rgcType,1}, spResponseCenter, spResponseSurround);
            
            %@JRG to move
            %             % For the EJ case.
            %
            %             if isa(ir.mosaic{rgcType},'rgcPhys')
            %                 rLinearSU = cell(length(ir.mosaic{rgcType}.sRFcenter));
            %                 for cellNum = 1:length(ir.mosaic{rgcType}.sRFcenter)
            %                     %                     rLinearSU{cellNum,1,1} = 3*responseLinear{cellNum,1}./max(responseLinear{cellNum,1}) + ir.mosaic{rgcType}.tonicDrive{cellNum,1};
            %                     rLinearSUTemp = 8.5*(responseLinear{cellNum,1} - ir.mosaic{rgcType}.tonicDrive{cellNum,1}) + ir.mosaic{rgcType}.tonicDrive{cellNum,1};
            %                     %                     rLinearSUTemp = 20*(responseLinear{cellNum,1} - ir.mosaic{rgcType}.tonicDrive{cellNum,1}) + ir.mosaic{rgcType}.tonicDrive{cellNum,1};
            %                     % NEED TO SUBSAMPLE TO GET CORRECT SPIKE RATE
            %                     rLinearSU{cellNum,1,1} = rLinearSUTemp;%(1:8:end);
            %                 end
            %                 clear responseLinear
            %                 responseLinear = rLinearSU;
            %             end
            
            % Store the linear response
            ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'responseLinear', responseLinear);
            
        end
    otherwise
        error('Unknown os type %s\n',osType);
end


