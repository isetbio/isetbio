function scene = sceneInterpolateW(scene, newWave, pLum)
% Wavelength interpolation for scene image data
%
% Syntax:
%	scene = sceneInterpolateW([scene], [newWave], [pLum ])
%
% Description:
%    Interpolate the wavelength dimension of a scene. By default, the
%    resampled scene has the same mean luminance as the original scene.
%
%    There are examples in the code. Type 'edit sceneInterpolateW' into
%    the Command Window to access.
%
% Inputs:
%    scene   - (Optional) Input scene structure. Default is to retrieve an
%              existing scene structure.
%    newWave - (Optional) Wavelengths of the output scene (interpolated).
%              Default is to query the user if not provided.
%    pLum    - (Optional) Boolean indicating whether or not to preserve the
%              luminance. Default is 1 (true).
%
% Outputs:
%    scene   - The modified scene structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - We do not think you should be extrapolating. We
%      therefore alert you if the newWave is outside of the range of the
%      current scene wavelengths.]
%    * N.B. The source contains executable examples of usage, which can be
%      accessed by typing 'edit sceneInterpolateW.m' in the command window.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/21/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match the Wiki

% Examples:
%{
    scene = sceneCreate;
    scene = sceneInterpolateW(scene, [400:10:700]);
%}
%{
    % Monochromatic scene
    scene = sceneCreate;
    scene = sceneInterpolateW(scene, 550);
    vcAddAndSelectObject(scene);
    sceneWindow;
%}
%{
    % Do not preserve luminance
    scene = sceneCreate;
    scene = sceneInterpolateW(scene, [400:2:700], 0);
    vcAddAndSelectObject(scene); sceneWindow;
%}

%% Initialize parameters
if notDefined('pLum'), pLum = 1; end
if notDefined('scene')
    scene = vcGetSelectedObject('scene');
elseif ~strcmp(sceneGet(scene, 'type'), 'scene')
    error('sceneInterpolationW: input structure not a scene');
end

% Check whether or not to show wait bar
showBar = ieSessionGet('wait bar');

% Note the current scene properties
curWave = sceneGet(scene, 'wave');

% If the user didn't send in new wavelengths, we ask.
if notDefined('newWave')
    handles = ieSessionGet('scene image handle');
    prompt = {'Start (nm)', 'Stop (nm)', 'Spacing (nm)'};
    def = {num2str(curWave(1)), num2str(curWave(end)), ...
             num2str(sceneGet(scene, 'binwidth'))};
    dlgTitle = 'Wavelength resampling';
    lineNo = 1;
    val = inputdlg(prompt, dlgTitle, lineNo, def);
    if isempty(val), return; end
    
    low = str2double(val{1});
    high = str2double(val{2});
    skip = str2double(val{3});
    
    if high >= low 
        newWave = low:skip:high;
    else
        ieInWindowMessage('Stop should be no larger than Start', ...
            handles, 5);
        return;
    end
end

% They are already the same.  Nothing to interpolate.  It happens.
if isequal(newWave(:),curWave(:)), return; end

waveSpectrum.wave = newWave;

if min(newWave) < min(curWave) ||  max(newWave) > max(curWave)
    error('Wavelength extrapolation not allowed.  Only interpolation');
end

%% Get current data and parameters
photons = sceneGet(scene, 'photons');
if ~isempty(photons) && pLum
    meanL = sceneGet(scene, 'mean luminance'); 
end

il = sceneGet(scene, 'illuminant');
if ~isempty(il)
    illuminantPhotons = sceneGet(scene, 'illuminant photons');
end

% Clear the current data before replacing.  This saves memory.
scene = sceneSet(scene, 'spectrum', waveSpectrum);
scene = sceneClearData(scene);

%% Interpolate the photons to the new wavelength sampling
if showBar, h = waitbar(0, 'Resampling wavelengths'); end
if ~isempty(photons)
    
    r = size(photons, 1);
    c = size(photons, 2);
    w = size(photons, 3);
    
    % Here is the extrapval problem
    if showBar, waitbar(0.3, h, 'Interpolating'); end
    
    % Not sure how big we should allow this to be
    if r * c * w < (640 * 640 * 31)
        % This is wavelength x space, rather than XW as usual
        [photons, row, col] = RGB2XWFormat(photons);     
        newPhotons = interp1(curWave, photons', newWave, 'linear')';
        newPhotons = XW2RGBFormat(newPhotons, row, col);
    else
        % Big data set condition, so we loop down the rows
        if ieSessionGet('gpu compute')
            newPhotons = zeros(...
                r, c, length(waveSpectrum.wave), 'gpuArray');
        else
            newPhotons = zeros(r, c, length(waveSpectrum.wave));
        end
        for rr = 1:r
            % Get a row
            pRow = squeeze(photons(rr, :, :))';
            % Interpolate all the columns in that row and put it in place
            newPhotons(rr, :, :) = interp1(curWave(:), pRow, ...
                waveSpectrum.wave(:), 'linear')';
        end
    end

    if showBar, waitbar(0.7, h, 'Compressing and storing'); end
    scene = sceneSet(scene, 'photons', newPhotons);
    
    % Calculate and store the scene luminance
    % scene = sceneSet(scene, 'luminance', sceneCalculateLuminance(scene));

end
if showBar, close(h); end

%% Now create the new illuminant. This over-writes the earlier one.
if ~isempty(il)
    
    % Interpolate the illuminant data
    newIlluminant = interp1(curWave, illuminantPhotons, newWave, 'linear');
    
    % Update the scene and illuminant spectrum.
    scene = sceneSet(scene, 'spectrum', waveSpectrum);
    scene = sceneSet(scene, 'illuminant spectrum', waveSpectrum);
    
    % Put in the new illuminant photons
    scene = sceneSet(scene, 'illuminant photons', newIlluminant);
end

%% Set and then adjust the luminance level
% For broadband scenes, we generally want to preserve the original mean
% luminance (stored in meanL) despite the resampling. In some cases, such
% as extracting a monochrome scene, we might not want to preserve the mean
% luminance.
if pLum && ~isempty(photons)
    scene = sceneAdjustLuminance(scene, meanL);
end

end