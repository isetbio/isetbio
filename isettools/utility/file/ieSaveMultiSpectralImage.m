function fullName = ieSaveMultiSpectralImage(fullName,mcCOEF,basis,comment,imgMean,illuminant,fov,dist)
%Save a Matlab data file containing a multi-spectral image.
%
%  fullName = ieSaveMultiSpectralImage(fullName,coef,basis,comment,imgMean,illuminant,fov,dist)
%
% Write the multispectral image file.  The variables are created using
% routines in the hypercube directory. 
% 
% Inputs:
%   mcCOEF  - coefficients (RGB format)
%   basis   - basis functions 
%   comment - 
%   imgMean - in some cases we remove the mean before creating the coeffs
%   illuminant structure
%     .wave  are wavelengths in nanometers
%     .data  are illuminant as a function of wavelength in energy units
%   fov     - Scene field of view
%   dist    - Distance to scene (should allow a depth map)
%
% fullName:  The full path to the output file is returned.
%
% The SPD of the data can be derived from the coefficients and basis
% functions using: 
%
%    spd = imageLinearTransform(mcCOEF,basis');
%
% See also: sceneToFile, hcBasis
%          
% Copyright ImagEval Consultants, LLC, 2005.

% TODO - Allow depth map for the scene, not just a distance

if notDefined('mcCOEF'),  error('Coefficients required');     end
if notDefined('basis'),   error('Basis function required.');  end
if notDefined('comment'), comment = sprintf('Date: %s\n',date); end %#ok<*NASGU>
if notDefined('illuminant'), error('Illuminant required'); end
if notDefined('fov'),     fov = 10; end    % 10 deg field of view is default
if notDefined('dist'),    dist = 1.2; end  % 1.2 meters distance is default
if notDefined('fullName')
    fullName = ...
        vcSelectDataFile('stayput','w','mat','Save multispectral data file.');
end
 
% Write out the matlab data file with the  information needed.
% Sometimes we save out data approximated using only the SVD
% Other times, we use a principal component method and have an image mean
% 12/2015: We added 'fov' and 'dist' as per NC.  
if notDefined('imgMean'), 
    save(fullName, 'mcCOEF', 'basis', 'comment', 'illuminant','fov','dist');
else
    save(fullName, 'mcCOEF', 'basis', 'imgMean', 'comment', 'illuminant','fov','dist');
end

end


