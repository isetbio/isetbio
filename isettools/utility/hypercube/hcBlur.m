function [hc, blur] = hcBlur(hc, sd)
% Blur each plane in a hypercube with conv2
%
% Syntax:
%   hc = hcBlur(hc, sd)
%
% Description:
%    Blur each plane in a hypercube with conv2. The spatial blur is a
%    gaussian created by fspecial with a size of sd.
%
% Inputs:
%    hc   - The hypercube
%    sd   - (Optional) The size of the symmetric Gaussian lowpass filter.
%           Default is 3.
%
% Outputs:
%    hc   - The hypercube
%    blur - The blur applied to each face of the hypercube
%

% History:
%    xx/xx/12       (c) Imageval 2012
%    12/05/17  jnm  Formatting

if notDefined('hc'), error('Hypercube data required.'); end
if notDefined('sd'), sd = 3; end

blur = fspecial('gaussian', [sd sd]);
nWave = size(hc, 3);
h = waitbar(0, 'Blurring');
for ii=1:nWave
    waitbar(ii / nWave, h);
    
    % Forces hypercube data to be doubles
    hc(:, :, ii) = conv2(double(hc(:, :, ii)), blur, 'same');
end
close(h);

end