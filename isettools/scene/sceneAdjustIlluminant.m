function scene = sceneAdjustIlluminant(scene, illEnergy, preserveMean)
% Adjust the current scene illuminant to the value in data
%
% Syntax:
%   scene = sceneAdjustIlluminant([scene], [illEnergy], [preserveMean])
%
% Description:
%    The scene radiance is scaled by dividing the current illuminant and
%    multiplying by the new illEnergy. The reflectance is preserved. If the
%    current scene has no defined illuminant, we assume that it has a D65
%    illumination. The scene luminance is preserved by this transformation.
%
%    There are examples in the code. Type 'edit sceneAdjustIlluminant' into
%    the Command Window to access.
% Inputs:
%    scene        - (Optional) A scene structure. Default is current scene.
%    illEnergy    - (Optional) A spectral data filename, or a vector
%                   matching the length of the scene wave. Defines the
%                   illuminant in energy units. Default is to query the
%                   user to select a file.
%    preserveMean - (Optional) A boolean indicating whether or not to
%                   preserve the mean illumination. Default is 'true'.
%
% Outputs:
%    scene        - The modified scene structure.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - Why do we have both warning and warndlg in the same
%      function? Which should be removed?]
%    * TODO: Determine illuminant estimation algorithm to insert.
%    * N.B. The source contains executable examples of usage, which can be
%      accessed by typing 'edit sceneAdjustIlluminant.m' in MATLAB's
%      command window.
%

% History:
%    xx/xx/10       Copyright ImagEval Consultants, LLC, 2010.
%    12/28/17  jnm  Formatting
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;  % Default is MCC under D65
    scene = sceneAdjustIlluminant(scene, 'Horizon_Gretag.mat');
    vcReplaceAndSelectObject(scene);
    sceneWindow;

    bb = blackbody(sceneGet(scene, 'wave'), 3000);
    scene = sceneAdjustIlluminant(scene, bb);
    vcReplaceAndSelectObject(scene);
    sceneWindow;

    bb = blackbody(sceneGet(scene, 'wave'), 6500, 'energy');
    figure;
    wave = sceneGet(scene, 'wavelength');
    plot(wave, bb)
    scene = sceneAdjustIlluminant(scene, bb);
    vcReplaceAndSelectObject(scene);
    sceneWindow;
%}

if notDefined('scene'), scene = vcGetObject('scene'); end
if notDefined('preserveMean'), preserveMean = true; end

% Make sure we have the illuminant data in the form of energy
wave = sceneGet(scene, 'wave');

if notDefined('illEnergy')
    % Read from a user-selected file
    fullName = vcSelectDataFile([]);
    illEnergy = ieReadSpectra(fullName, wave);
elseif ischar(illEnergy)
    % Read from the filename sent by the user
    fullName = illEnergy;
    if ~exist(fullName, 'file')
        error('No file %s\n', fullName);
    else
        illEnergy = ieReadSpectra(fullName, wave);
    end
else
    % User sent numbers. We check for numerical validity next.
    fullName = '';
end

% We check the illuminant energy values.
if max(illEnergy) > 10 ^ 5
    % Energy is not this big.
    warning(['Illuminant energy values are high; may be photons, ' ...
        'not energy.'])
elseif isequal(max(isnan(illEnergy(:))), 1) ...
        || isequal(min(illEnergy(:)), 0)
    warndlg(['NaNs or zeros present in proposed illuminant over this ' ...
        'wavelength range. No transformation applied.']);
    pause(3);
    return;
end

% Start the conversion
curIll = sceneGet(scene, 'illuminant photons');
if isempty(curIll)
    % [Note: XXX - We treat this as an opportunity to create an illuminant,
    % as in sceneFromFile (or vcReadImage). Assume the illuminant is D65.
    % Lord knows why. Maybe we should do an illuminant estimation algorithm
    % here. -- TODO: Determine illuminant estimation algorithm to insert.]
    disp('Old scene. Creating D65 illuminant')
    wave = sceneGet(scene, 'wave');
    curIll = ieReadSpectra('D65', wave);  % D65 in energy units
    scene = sceneSet(scene, 'illuminant energy', curIll);
    curIll = sceneGet(scene, 'illuminant photons');
end

% Current mean luminance may be preserved
mLum = sceneGet(scene, 'mean luminance');
if isnan(mLum) 
    [lum, mLum] = sceneCalculateLuminance(scene);
    scene = sceneSet(scene, 'luminance', lum);
end

% Converts illEnergy to illPhotons. Deals with different illuminant
% formats. If preserve reflectance or not, do slightly different things.
curIll = double(curIll);
switch sceneGet(scene, 'illuminant format')
    case 'spectral'
        % In this case the illuminant is a vector. We convert to photons
        illPhotons = Energy2Quanta(illEnergy, wave);

        % Find the multiplier ratio 
        illFactor = illPhotons ./ curIll;
        
        % Adjust the radiance data and the illuminant by the illFactor
        % This preserves the reflectance.
        % [Note: XXX - Don't skip changing the illuminant (do change it!)]
        skipIlluminant = 0;
        scene = sceneSPDScale(scene, illFactor, '*', skipIlluminant);
        
    case 'spatial spectral'
        if ~isequal(size(illEnergy), size(illEnergy))
            error('Spatial spectral illuminant size mis-match');
        end
        
        [newIll, r, c] = RGB2XWFormat(illEnergy);
        newIll = Energy2Quanta(wave, newIll');
        newIll = XW2RGBFormat(newIll', r, c);
        
        % Get the scene radiance
        photons = sceneGet(scene, 'photons');
        
        % Divide the radiance photons by the current illuminant and then
        % multiply by the new illuminant. These are the radiance photons
        % under the new illuminant. This preserves the reflectance.
        photons = (photons ./ curIll) .* newIll;
        
        % Set the new radiance back into the scene
        scene = sceneSet(scene, 'photons', photons);
        
        % Set the new illuminant back into the scene
        scene = sceneSet(scene, 'illuminant photons', newIll);
end

% Make sure the mean luminance is unchanged. (preseveMean defaults to true)
if preserveMean, scene = sceneAdjustLuminance(scene, mLum); end

scene = sceneSet(scene, 'illuminant comment', fullName);

return;
