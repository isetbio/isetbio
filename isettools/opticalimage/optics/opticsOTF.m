function oi = opticsOTF(oi, scene)
% Apply the opticalImage OTF to the photon data
%
% Syntax:
%   oi = opticsOTF(oi, scene);
%
% Description:
%    The optical transform function (OTF) associated with the optics in the
%    OI is calculated and applied to the scene data. This function is
%    called for shift-invariant and diffraction-limited models. It is not
%    called for the ray trace approach.
%
%    The spatial (frequency) support of the OTF is computed from the OI
%    information.
%
%    The OTF data are not stored or returned. The OTF can be quite large.
%    It represents every spatial frequency in every waveband. So we compute
%    the OTF and apply it on the fly, without ever representing the whole
%    OTF (support and wavelength).
%
%    The programming issues concerning using Matlab to apply the OTF to the
%    image (rather than convolution in the space domain) are explained both
%    here and in the script s_FFTinMatlab.
%
%    In the future, we may permit saving the OTF in the OI structure by
%    setting the otfSaveFlag to true (1). At present, that flag is not
%    implemented. But computers seem to be getting bigger and faster.
%
% Inputs:
%    oi    - Struct. An optical image structure.
%    scene - (Optional) Struct. A scene structure. Default uses vcGetObject
%            to retrieve an existing scene structure.
%
% Outputs:
%    oi    - Struct. The modified optical image structure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * TODO: Implement otfSaveFlag
%
% See Also:
%   s_FFTinMatlab, oiCalculateOTF, oiCompute
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    03/09/18  jnm  Formatting
%    04/07/18  dhb  Fixe opticsOTF example so it runs.
%                   Remove examples from helper functions in this file.
%                   They can't work because you can't call those helper
%                   functions from outside of the top level function.
%    10/04/18  npc  Handle custom oi padding
%    06/27/19  JNM  Minor formatting adjustments

% Examples:
%{
   oi = oiCreate('human');
   oi = opticsOTF(oi);
%}

if notDefined('oi'), error('Optical image required.'); end
if notDefined('scene'), scene = vcGetObject('scene'); end

opticsModel = oiGet(oi, 'optics model');

switch lower(opticsModel)
    case {'skip', 'skipotf'}
        return
    case {'dlmtf', 'diffractionlimited'}
        oi = oiApplyOTF(oi, scene);
    case {'shiftinvariant', 'custom', 'humanotf'}
        oi = oiApplyOTF(oi, scene, 'mm');
    otherwise
        error('Unknown OTF method');
end

end

%-------------------------------------------
function oi = oiApplyOTF(oi, scene, unit)
% Calculate and apply the otf waveband by waveband
%
% Syntax:
%   oi = oiApplyOTF(oi, method, unit);
%
% Description:
%    We calculate the OTF every time, never saving it, because it can take
%    up a lot of space and is not that hard to calculate. Also, any change
%    to the optics properties would make us recompute the OTF, and keeping
%    things synchronized can be error prone.
%
% Inputs:
%    oi    - Struct. An optical image structure.
%    scene - Struct. A scene structure.
%    unit  - (Optional) String. A string describing what units to use.
%            Default is 'cyclesPerDegree'.
%
% Outputs:
%    oi    - Struct. The modified optical image structure.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/09/18  jnm  Formatting

if notDefined('unit'), unit = 'cyclesPerDegree'; end
wave = oiGet(oi, 'wave');

% Pad the optical image to allow for light spread. Also, make sure the row
% and col values are even.

%% Determine padSize, padValue and rectRadius
% A non-empty returned rectRadius indicates that the user-supplied
% padSizeDegs will result in a pad size that is smaller than the default,
% which is 1/8 of the oi width. In this case, we pad using the default pad
% size and call oiCrop using the rectRadius value after convolution with
% the PSF. This is done at the end of this function.
[padSize, padValue, rectRadius] = oiPadParams(oi);

sDist = sceneGet(scene, 'distance');
oi = oiPad(oi, padSize, padValue, sDist);

% See s_FFTinMatlab to understand the logic of the operations here. We used
% to do this one wavelength at a time. But this could cause dynamic range
% problems for ieCompressData. So, for now we are experimenting with
% filtering one at a time but stuffing the whole data set in at once.

% Get the current data set. It has the right size. We over-write it below.
p = oiGet(oi, 'photons');
otfM = oiCalculateOTF(oi, wave, unit);

for ii = 1:length(wave)
    % img = oiGet(oi, 'photons', wave(ii));
    img = p(:, :, ii);
    % figure(1);
    % imagesc(img);
    % colormap(gray);

    % For diffraction limited we calculate the OTF. For other optics models
    % we look up the stored OTF. Remember, DC is in the (1, 1) position.
    otf = otfM(:, :, ii);
    % vcNewGraphWin;
    % mesh(fftshift(otf));
    % otf(1, 1)

    % Put the image center in (1, 1) and take the transform.
    imgFFT = fft2(fftshift(img));
    % figure(1);
    % imagesc(abs(imgFFT));
    % figure(2);
    % imagesc(abs(otf));
    % colormap(gray)

    % Multiply the transformed otf and the image.
    % Then invert and put the image center in  the center of the matrix
    filteredIMG = abs(ifftshift(ifft2(otf .* imgFFT)));

    % Sometimes we had annoying complex values left after this filtering.
    % We got rid of it by an abs() operator. It should never be there. But
    % we think it arises because of rounding error. We haven't seen this in
    % years, however.
    % figure(1);
    % imagesc(abs(filteredIMG));
    % colormap(gray)
    p(:, :, ii) = filteredIMG;
end

% Put all the photons in at once.
oi = oiSet(oi, 'photons', p);

% Non-empty returned rectRadius. Crop computed oi accordingly.
if (~isempty(rectRadius))
    s = oiGet(oi, 'size');
    oiCenter = round(s / 2);
    cropRect = [oiCenter(1) - rectRadius, oiCenter(2) - rectRadius, ...
        rectRadius * 2, rectRadius * 2];
    oi = oiCrop(oi,cropRect);
end

end
