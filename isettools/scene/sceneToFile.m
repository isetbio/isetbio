function varExplained = sceneToFile(fname, scene, bType, mType)
% Write scene data in the hyperspectral and multispectralfile format
%
% Syntax:
%	varExplained = sceneToFile(fname, scene, [bType], [mType])
%
% Description:
%    If the cFlag is empty, it saves a file containing photons, wave, 
%    illuminant structure and a comment.
%
%    If the cFlag is a value (double), the function builds a linear model
%    basis to represent (and compress) the photon data. It saves the linear
%    model, model coefficients, illuminant structure, and a comment. The
%    linear model format removes the mean of the photons, builds a linear
%    model, and stores the mean, linear model, and coefficients for each
%    pixel.
%
% Inputs:
%    fname        - The full name of the output file
%    scene        - ISETBIO scene structure 
%	 bType        - (Optional) Basis calculation type. Default 0.99
%                   Empty - No compression, just save photons, wave,
%                            comment, illuminant
%                   0 - 1 - A value between 0 and 1 will specify the
%                            fraction of variance explained by the linear
%                            model compresssion (Default 0.99)
%                   > 1   - Integer specifying a number of basis functions
%	 mType        - (Optional) Mean computation. Default 'canonical'.
%                   Remember to remove the mean 'meansvd' or not
%                   'canonical' before calculating svd.
% Outputs:
%	 varExplained - Fraction of variance explained by the linear model
%
% Optional key/value pairs:
%    None.
%
% Examples of usage are provided in the source code, and can be
% accessed by typing 'edit sceneToFile.m' in MATLAB's command window.
%
% Notes:
%    * [NOTE: XXX - Add depth image as potential output, not just dist.]
%    * N.B. The source contains executable examples of usage, accessible by
%      typing 'edit sceneToFile.m' in the command window.
%

% History:
%    xx/xx/13       (c) Imageval Consulting, LLC 2013
%    12/20/17  jnm  Formatting & fix default assignments
%    01/25/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    scene = sceneCreate;
    vcAddAndSelectObject(scene);
    sceneWindow;
    fname = fullfile(tempdir,'deleteMePlease');
    sceneToFile(fname, scene, 0.999);
    scene2 = sceneFromFile(fname, 'multispectral');
    ieAddObject(scene2);
    sceneWindow;
    sceneToFile(fname, scene, []);
%}

if notDefined('fname'), error('Need output file name for now'); end
if notDefined('scene'), error('scene structure required'); end
if notDefined('bType'), bType = 0.99; end % See hcBasis
if notDefined('mType'), mType = 'canonical'; end % See hcBasis

% We need to save the key variables
photons = sceneGet(scene, 'photons');
wave = sceneGet(scene, 'wave');
illuminant = sceneGet(scene, 'illuminant');
fov = sceneGet(scene, 'fov');
dist = sceneGet(scene, 'distance');
comment = sprintf('Scene: %s', sceneGet(scene, 'name'));

if isempty(bType)
    % No compression.
    save(fname, 'photons', 'wave', 'comment', 'illuminant', 'fov', 'dist');
    varExplained = 1;
else
    % Figure out the basis functions using hypercube computation
    photons = photons(1:3:end, 1:3:end, :);
    [imgMean, basisData, ~, varExplained] = hcBasis(photons, bType, mType);
    clear photons;
    
    % Plot the basis functions
    %   wList = sceneGet(scene, 'wave');
    %   vcNewGraphWin;
    %   for ii = 1:size(basisData, 2)
    %       plot(wList, basisData(:, ii));
    %       hold on
    %   end   
    
    photons = sceneGet(scene, 'photons');
    [photons, row, col] = RGB2XWFormat(photons);
    switch mType
        case 'canonical'
            coef = photons * basisData;
            
        case 'meansvd'
            photons = photons - repmat(imgMean, row * col, 1);
            coef = photons * basisData;
            
        otherwise
            error('Unknown mType: %s\n', mType);
    end
    coef = XW2RGBFormat(coef, row, col);

    % Save the coefficients and basis
    basis.basis = basisData;
    basis.wave = wave;
    ieSaveMultiSpectralImage(fname, coef, basis, comment, imgMean, ...
        illuminant, fov, dist);

end

end  % End function
