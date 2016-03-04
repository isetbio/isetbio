function ir = irCompute(ir, outerSegment, varargin)
% Computes the rgc mosaic responses to an arbitrary stimulus.
%
% The responses for all the mosaics attached to the ir object.
%
% For each mosaic, the linear response is computed by spatial convolution
% of the center and surround RFs. Then, the temporal responses for the
% center and surround is computed, separably. This stage of the computation
% is stored in responseLinear.  
%
% The linear response is then sent to the irComputeSpikes routine.  Each
% rgc type will have its own compute method.  For all types the initial
% linear method is carried out in here.
%
% Inputs: 
%  ir: inner retina object, 
%  outersegment: ojbect
%
% Outputs: 
%  ir: the inner retina object with responses attached 
%
% Example:
%   ir.compute(identityOS);
%   irCompute(ir,identityOS);
%
% See also: rgcGLM, irComputeSpikes, 
%
% JRG (c) isetbio team

%% PROGRAMMING TODO
% We plan to move the linear computation to a separate routine,
% irComputeContinuous (or Linear) and carry out the first part of this code
% in that routine, rather than as the switch statement that is here now
%

%% Parse inputs
p = inputParser;
p.CaseSensitive = false;

p.addRequired('ir',@(x) ~isempty(validatestring(class(x),{'ir','irPhys'})));
p.addRequired('outerSegment',@(x) ~isempty(validatestring(class(x),{'osIdentity','osLinear','osBioPhys'})));

p.parse(ir,outerSegment,varargin{:});

%% Get the input data

% Possible osTypes are osIdentity, osLinear, and osBiophys
% Only osIdentity is implemented now.
osType = class(outerSegment);
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
        
    case {'osLinear'}
        %% Linear OS
        error('Not yet implemented');
    case {'osBioPhys'}
        %% Full biophysical os
        error('Not yet implemented');
    otherwise
        error('Unknown os type %s\n',osType);
end


%% Linear computation

% This will end up begin a call toe irComputeContinuous

% Determine the range of the rgb input data
spTempStim = osGet(outerSegment, 'rgbData');
range = max(spTempStim(:)) - min(spTempStim(:));

% Special case. If the ir class is rgcPhys, which we use for
% validation. But in general, this is not the case. There is a
% scientific question about this 10.  We need JRG and EJ to resolve
% the spiking rate.
if isequal(class(ir),'irPhys'),   spTempStim = spTempStim./range;
else                    spTempStim = 10*spTempStim./range;
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
    
    % Convolve with the temporal impulse response
    responseLinear = ...
        fullConvolve(ir.mosaic{rgcType,1}, spResponseCenter, spResponseSurround);
    
    % Store the linear response
    ir.mosaic{rgcType} = mosaicSet(ir.mosaic{rgcType},'responseLinear', responseLinear);
    
end

% Compute spikes for each trial
switch class(ir.mosaic{rgcType})
    case {'rgcLinear','rgcPhys'};
        % No nonlinear response
        error('Not yet implemented');s
    otherwise
        nTrials = ir.mosaic{rgcType}.numberTrials;
        for itrial = 1:nTrials
            ir = irComputeSpikes(ir);
        end
end


