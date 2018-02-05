function fullName = ieSaveMultiSpectralImage(fullName, mcCOEF, basis, ...
    comment, imgMean, illuminant, fov, dist)
% Save a Matlab data file containing a multi-spectral image.
%
% Syntax:
%   fullName = ieSaveMultiSpectralImage([fullName], mcCOEF, basis, ...
%       [comment], [imgMean], illuminant, [fov], [dist])
%
% Description:
%    Write the multispectral image file. The variables are created using
%    routines in the hypercube directory. 
%
%    Examples are located within the code. To access the examples, type
%    'edit ieSaveMultiSpectralImage.m' into the Command Window.
% 
% Inputs:
%    fullName   - (Optional) full file name and path. Default is looked up.
%    mcCOEF     - coefficients (RGB format)
%    basis      - basis functions 
%    comment    - (Optional) If comment is not provided, use current date
%    imgMean    - (Optional) in some cases we remove the mean before 
%                 creating the coeffsilluminant structure. Default is not
%                 providing an imgMean and not passing it along.
%        .wave  are wavelengths in nanometers
%        .data  are illuminant as a function of wavelength in energy units
%    illuminant - The illuminant
%    fov        -(Optional) Scene field of view. Default is 10 degrees
%    dist       - (Optional) Distance to scene (should allow a depth map).
%                 Default is 1.2 meters.
%
% Outputs:
%    fullName - The full path to the output file is returned.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * The SPD of the data can be derived from the coefficients and basis
%      functions using: 
%          spd = imageLinearTransform(mcCOEF, basis');
%    * TODO: Allow depth map for the scene, not just a distance
%
% See Also:
%    sceneToFile, hcBasis
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/27/17  jnm  Formatting
%    11/29/17  jnm  Add note & Example
%    01/26/18  jnm  Formatting update to match the Wiki.

% Examples:
%{
    % Read in the scene
    fName = fullfile(isetbioDataPath, 'images', 'multispectral', ...
        'StuffedAnimals_tungsten-hdrs');
    scene = sceneFromFile(fName, 'multispectral');

    % Compress the hypercube requiring only 95% of the var explained
    vExplained = 0.95;
    [imgMean, imgBasis, coef] = hcBasis(sceneGet(scene, 'photons'), ...
        vExplained);

    % Save the data 
    wave = sceneGet(scene, 'wave');
    basis.basis = imgBasis;
    basis.wave = wave;

    comment = 'Compressed using hcBasis with imgMean)';

    illuminant = sceneGet(scene, 'illuminant');
    % illuminant.wavelength = scene.spectrum.wave;
    % illuminant.data = scene.illuminant.data;
    oFile = fullfile(isetRootPath, 'deleteMe.mat');
    ieSaveMultiSpectralImage(oFile, coef, basis, comment, imgMean, ...
        illuminant);
    delete(oFile);
%}

if notDefined('mcCOEF'), error('Coefficients required');     end
if notDefined('basis'), error('Basis function required.');  end
if notDefined('comment'), comment = sprintf('Date: %s\n', date); end
if notDefined('illuminant'), error('Illuminant required'); end
if notDefined('fov'), fov = 10; end     % 10 deg field of view is default
if notDefined('dist'), dist = 1.2; end  % 1.2 meters distance is default
if notDefined('fullName')
    fullName = ...
        vcSelectDataFile('stayput', 'w', 'mat', ...
            'Save multispectral data file.');
end 
if notDefined('imgMean')
    save(fullName, 'mcCOEF', 'basis', 'comment', 'illuminant', 'fov', ...
        'dist');
else
    save(fullName, 'mcCOEF', 'basis', 'imgMean', 'comment', ...
        'illuminant', 'fov', 'dist');
end

end
