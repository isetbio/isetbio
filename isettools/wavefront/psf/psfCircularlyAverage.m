function outPSF = psfCircularlyAverage(inPSF)
% Circularly average the provided PSF
%
% Syntax:
%   outPSF = psfCircularlyAverage(inPSF)
%
% Description:
%    As the name suggests. The PSF output volume is scaled to match the PSF
%    input volume.
% 
% Inputs:
%    inPSF  - Point-Spread Function
%
% Outputs:
%    outPSF - Circularly-averaged Point-Spread Function
%

% History:
%    07/19/07  dhb  Wrote it.
%    12/22/09  dhb  Fix bug in how peakRow and peakCol are computed.
%    12/22/09  dhb  Make computation a little more fine grained.
%    07/23/12  dhb  Match out volume to in volume.
%    11/13/17  jnm  Comments & Formatting
%

% Examples:
%{
    % This example only shows the circularly averaged PSF for L-cones
    theZernikeCoeffs = importdata('autrusseauStandardObserver.txt');
    % Cone sensitivities and equal energy weighting spectrum
    load('T_cones_ss2');
    conePsfInfo.S = S_cones_ss2;
    conePsfInfo.T = T_cones_ss2;
    conePsfInfo.spdWeighting = ones(conePsfInfo.S(3),1);

    wls = SToWls([400 10 31]);
    wvf0 = wvfCreate;

    % Set important parameters - Autrusseau std. observer
    wvf0 = wvfSet(wvf0,'measured pupil size',6);
    wvf0 = wvfSet(wvf0,'calc pupil size',6);
    wvf0 = wvfSet(wvf0,'zcoeffs',theZernikeCoeffs(:,1));
    wvf0 = wvfSet(wvf0,'measured wavelength',570);
    wvf0 = wvfSet(wvf0,'calc wavelengths',wls);
    wvf0 = wvfSet(wvf0,'calc cone psf info',conePsfInfo);
    wvf0 = wvfSet(wvf0,'number spatial samples',497);
    sce = sceCreate(wls,'none');
    wvf0 = wvfSet(wvf0,'sce params',sce);

    wvfParams1 = wvf0;
    wvfParams1 = wvfComputePSF(wvfParams1);
    conePsf1 = wvfGet(wvfParams1,'cone psf');

    lpsf = conePsf1(:,:,1);
    lpsf = psfCircularlyAverage(lpsf);
%}

% Define quantization. Four was used in early code, but 1 makes more sense.
quantizationFactor = 1;

% Make a circularly symmetric version of average optics.
[m, n] = size(inPSF);
if (n ~= m)
    error('Input must be a square matrix');
end
nLinearPixels = m;

[peakRow, peakCol] = psfFindPeak(inPSF);
radiusMat = MakeRadiusMat(nLinearPixels, nLinearPixels, peakCol, peakRow);
outPSF = zeros(nLinearPixels, nLinearPixels);
nBands = round(nLinearPixels / quantizationFactor);
radii = linspace(0, 0.75 * nLinearPixels, nBands);
for q = 1:length(radii) - 1
    index = find(radiusMat >= radii(q) & radiusMat < radii(q + 1));
    if (~isempty(index))
        outPSF(index) = mean(inPSF(index));
    end
end
outPSF = sum(inPSF(:)) * outPSF / sum(outPSF(:));
