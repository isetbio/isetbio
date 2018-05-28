function scene = sceneFromBasis(sceneS)
% Create a scene from a structure with basis functions (linear model)
%
% Syntax:
%	scene = sceneFromBasis(sceneS)
%
% Description:
%    Using a linear model, create a scene from a structure (sceneS) with
%    basis functions.
%
%    This function is used with the remote data toolbox, which returns a
%    Matlab scene formatted as a set of basis functions and coefficients.
%
%    There are examples contained in the code. To access, simply type 'edit
%    sceneFromBasis.m' into the Command Window.
%
% Inputs:
%   sceneS - a structure that has mcCOEF, basis, and illuminant
%
% Output:
%   scene  - A scene structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Unsure if example is functioning as desired.]
%    * N.B. The source contains executable examples of usage, which can be
%      accessed by typing 'edit sceneFromBasis.m' in the command window.
%

% History:
%    xx/xx/xx  BW   ISETBIO Team
%    12/22/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    % ETTBSkip - Error with gradleFetchArtifact?
    rd = RdtClient('isetbio');
    rd.crp('/resources/scenes/multiband/scien/2004');
    data = rd.readArtifact('AsianFemale_2', 'type', 'mat');
    scene = sceneFromBasis(data);
    ieAddObject(scene);
    sceneWindow;
%}

%%
if ~isfield(sceneS, 'mcCOEF'), error('mcCOEF required'); end
if ~isfield(sceneS, 'basis'), error('basis required'); end
imgMean = [];
if isfield(sceneS, 'imgMean'), imgMean = single(sceneS.imgMean); end

%% The image data should be in units of photons
photons = imageLinearTransform(sceneS.mcCOEF, sceneS.basis.basis');

if ~isempty(imgMean)
    % The saved function was calculated using principal components, not
    % just the SVD. Hence, the mean is stored and we must add it into the
    % computed image.
    [photons, r, c] = RGB2XWFormat(photons);
    photons = repmat(imgMean(:), 1, r * c) + photons';
    photons = XW2RGBFormat(photons', r, c);
end

% Force photons to be positive. Sometimes there is a tiny negative value
photons = max(photons, 0);
% vcNewGraphWin; imageSPD(photons, basis.wave);

% These lines are left in because there must be different file types out
% there somewhere. Sometimes we stored the mean, and sometimes we didn't.
% But apparently it is rare.

scene = sceneCreate('multispectral');
scene = sceneSet(scene, 'wavelength', sceneS.basis.wave);
scene = sceneSet(scene, 'photons', photons);

% Deal with the illuminant
illuminant = illuminantModernize(sceneS.illuminant);
scene = sceneSet(scene, 'illuminant', illuminant);

end
