function [hc, blur] = hcBlur(hc, sd)
% Blur each plane in a hypercube with conv2
%
% Syntax:
%   [hc, blur] = hcBlur(hc, sd)
%
% Description:
%    Blur each plane in a hypercube with conv2. The spatial blur is a
%    gaussian created by fspecial with a size of sd.
%
%    Examples are located within the code. To access the examples, type
%    'edit hcBlur.m' into the Command Window.
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
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/12       (c) Imageval 2012
%    12/05/17  jnm  Formatting
%    12/27/17   BW  Fixed blur bug and changed return
%    01/23/18  dhb  Fix example so it runs.
%    01/26/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    mcc = sceneCreate;
    hc = sceneGet(mcc, 'photons');
    [hc, blurFunction] = hcBlur(hc, 5);
    imageSPD(hc,sceneGet(mcc, 'wave'));
    vcNewGraphWin;
    mesh(blurFunction);
    colormap(jet)
%}

if notDefined('hc'), error('Hypercube data required.'); end
if notDefined('sd'), sd = 3; end

blur = fspecial('gaussian', 3 * sd, sd);
nWave = size(hc, 3);
h = waitbar(0, 'Blurring');
for ii = 1:nWave
    waitbar(ii / nWave, h);
    
    % Forces hypercube data to be doubles
    hc(:, :, ii) = conv2(double(hc(:, :, ii)), blur, 'same');
end
close(h);

end