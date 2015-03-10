function [quantImg,quantizationError] = analog2digital(ISA,method)
%Analog-to-digital conversion for voltage data in a sensory array
%  
%  [qImage,qError] = analog2digital(ISA,method)
%
%  Various quantization schemes are implemented. This routine calculates
%  the quantized image and the quantization error, if requested.. 
%   
% Example:
%
% [qImage,qError] = analog2digital(ISA,method);
%
% Copyright ImagEval Consultants, LLC, 2005.

%%
if notDefined('ISA'),    ISA = vcGetObject('sensor'); end
if notDefined('method'), method = sensorGet(ISA,'quantizationMethod'); end

%% Get voltage data and range
voltageSwing = pixelGet(ISA.pixel,'voltageSwing');
img          = sensorGet(ISA,'volts'); 
if isempty(img), error('No voltage image'); end

%% Apply method
switch lower(method)
    case {'analog'}
        quantImg = img;
        if nargout == 2
            quantizationError = zeros(size(img)); 	% [mV]
        end
        
    case {'lin', 'linear'}
        nBits = sensorGet(ISA,'nbits'); 
        if isempty(nBits)
            nBits = 8; 
            warning('ISET:Quantization0','Assuming %d bits.',nBits);
        end
        quantizationStep = voltageSwing / (2^nBits);	    % [mV/DN]
        quantImg = round(img/quantizationStep);             % [DV]
        if nargout == 2
            quantizationError = img - (quantImg * quantizationStep); 	% [mV]
        end
    otherwise
        warning('ISET:Quantization1','Unknown quantization method.')
end

end