function bMatrix = colorBlockMatrix(wList, extrapVal)
% Create a matrix to render spd data into RGB
%
% Syntax:
%   bMatrix = colorBlockMatrix(wList, extrapVal)
%
% Description:
%    We render spectral data in the scene and optical image windows as RGB
%    images. The matrix returned by this routine is used to calculate R, G
%    and B values from the SPD. The columns of the returned matrix define
%    how to sum across the wavebands.
%
%    There are many ways one could do an spd -> RGB rendering. This one is
%    simple, and uses a linear mapping that is implemented by the matrix
%    returned by this routine. By default, the wavelengths from 400-490 add
%    to the blue channel, from 500-570 add to the green channel, and
%    580-700 add to the red channel. Wavelengths outside of this band, by
%    default, do not contribute.
%
%    It is possible to generate a contribution from outside the band, say
%    in the infrared. When we are trying to visualize IR, it is useful to
%    set a value of 0.1 or 0.2.
%
%    This function contains examples of usage inline. To access these, type
%    'edit colorBlockMatrix.m' into the Command Window.
%
% Inputs:
%    wList     - Vector. The list of wavelengths in the SPD to be rendered.
%    extrapval - (Optional) Numeric. The amount contributed outside the
%                visible band. Default is 0.
%
% Outputs:
%    bMatrix   - Matrix. The color matrix used for rendering a photon
%                spectrum, where photonSPD is a column vector, following
%                the following format:
%                   displayRGB = photonSPD * bMatrix
%
% Optional key/value pairs:
%    whiteSPD  - Vector. By default, an equal photon spectrum is rendered
%                as white (1, 1, 1). To render another photon spectrum as
%                (1, 1, 1), send in this vector.
%
% Notes:
%    * [NOTE - XXX: We used to set this with ieSessionSet and manage it
%      with ieSessionGet and matlab setpref/getpref. Now, not so much. This
%      code is left here as a reminder that we might reconsider. (WP option
%      section, currently hardcoded to d65).]
%    * [Note - XXX: wp switch section - We should think about whether we
%      want whiteSPD to be persistent]
%    * [Note - XXX: wp switch section - We used to set this with
%      ieSessionSet and manage it with ieSessionGet and matlab
%      setpref/getpref. Now, not so much. This code is left here as a
%      reminder that we might reconsider.]
%    * TODO - If wList a single number, then we interpret it as one of the
%      common visible wavelength ranges. This feature is here for backwards
%      compatibility. But I am trying to eliminate all such calls. -
%      Eliminate calls of this nature.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC.
%    10/30/17  jnm  Comments & formatting
%    07/11/19  JNM  Formatting update

% Examples:
%{
    wList = [400:5:700];
    bMatrix = colorBlockMatrix(wList);
    figure;
    plot(wList, bMatrix)
%}
%{
    wList = [400:5:900];
    bMatrix = colorBlockMatrix(wList, 0.1);
    figure;
    plot(wList, bMatrix)
%}
%{
    % The spectrum that will be rendered as white (1, 1, 1) is equal photon
    % This can be calculated as:
    wList = [400:5:700];
    bMatrix = colorBlockMatrix(wList);
    whiteSPDPhotons = [1, 1, 1] * pinv(bMatrix);
    figure;
    plot(wList, whiteSPDPhotons)
%}

if notDefined('wList')
    error('You must specify wavelengths');
end
if notDefined('extrapVal')
    extrapVal = 0.0;
end

% If wList a single number, then we interpret it as one of the common
% visible wavelength ranges. This feature is here for backwards
% compatibility. But I am trying to eliminate all such calls
if length(wList) == 1
    % warning('Old wList format.')
    if wList == 31
        wList = 400:10:700;
    elseif wList == 371
        wList = 370:730;
    elseif wList == 37
        wList = 370:10:730;
    end
end

% The default block matrix function is defined over 400:10:700 range
defaultW = (400:10:700);

% Equal energy looks pink this way. Maybe we should adjust?
% Equal photons looks white this way.
% D65 looks kind of OK.
b = 10;
g = 8;
r = 31 - b - g;
defaultMatrix = [zeros(1, b), zeros(1, g),  ones(1, r); ...
    zeros(1, b),  ones(1, g), zeros(1, r); ...
    ones(1, b), zeros(1, g), zeros(1, r)]';

% Adjust for any differences in the wave list
if isequal(wList(:), (400:10:700)')
    % Set the default matrix columns sum to 1
    d = sum(defaultMatrix);
    bMatrix = defaultMatrix * diag(1 ./ d);
else
    % Adjust the matrix to match the default over 400-700 but be a small
    % value in the infrared. The default is 0. But it could be a non-zero
    % value the user sends in.
    bMatrix = zeros(length(wList), 3);
    for ii = 1:3
        bMatrix(:, ii) = interp1(defaultW(:), defaultMatrix(:, ii), ...
            wList(:), 'linear', extrapVal);
    end
    d = sum(bMatrix);
    bMatrix = bMatrix * diag(1 ./ d);
end

% As currently configured, the default sets equal photon input to (1, 1, 1)
% That is
%   ones(1, length(wList)) * bMatrix = (1, 1, 1)
% Examples below show how to change from default in which equal photon
% equal to (1, 1, 1)
% We should think about whether we want whiteSPD to be persistent
% We used to set this with ieSessionSet and manage it with ieSessionGet and
% matlab setpref/getpref. Now, not so much. This code is left here as a
% reminder that we might reconsider.
wp = 'd65';
switch lower(wp)
    case 'ee'
        % Make equal energy (1, 1, 1)
        ee = ones(length(wList), 1);
        eePhotons = Energy2Quanta(wList, ee);
        whiteSPD  = eePhotons / max(eePhotons);
    case 'd65'
        % Make D65 (1, 1, 1)
        d65Photons = blackbody(wList, 6500, 'photons');
        whiteSPD = d65Photons / max(d65Photons);
    otherwise
        % ep (equal photon spd) goes here
        % Default is equal photon count maps to (1, 1, 1)
        whiteSPD = [];
end

% Premultiply by the whiteSPD, if it was sent in. In this case, an SPD
% whose photons are given by whiteSPD will be mapped to (1, 1, 1).
if ~notDefined('whiteSPD')
    bMatrix = diag(1 ./ whiteSPD(:)) * bMatrix;
end

end