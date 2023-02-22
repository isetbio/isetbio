function [impulseResponse, temporalSupport] = impulseResponseFromOneSidedTransferFunction(oneSidedTransferFunction, tfSupport)

    nSamples = length(oneSidedTransferFunction);
    oneSidedTransferFunction(2:end-1) = oneSidedTransferFunction(2:end-1)/2;    % Divide by 2 to correct for amplitude

    % Add second size (flipped and conjugated)
    doubleSidedTransferFunction = [oneSidedTransferFunction fliplr(conj(oneSidedTransferFunction(2:nSamples)))];
    impulseResponse = ifft(doubleSidedTransferFunction, 'symmetric');

    df = (tfSupport(2)-tfSupport(1));
    nF = numel(tfSupport);
    dt = 0.5/(df*nF);
    temporalSupport = (0:(numel(impulseResponse)-1))*dt;
end