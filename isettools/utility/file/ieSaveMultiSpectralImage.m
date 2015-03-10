function fullName = ieSaveMultiSpectralImage(fullName,mcCOEF,basis,comment,imgMean,illuminant)
%Save a Matlab data file containing a multi-spectral image.
%
%  fullName = ieSaveMultiSpectralImage(fullName,coef,basis,comment,imgMean,illuminant)
%
% The image is created using routines in the multicapture directory.
% 
% Input arguments
%   mcCOEF  - coefficients (RGB format)
%   basis   - basis functions 
%   comment - 
%   imgMean - in some cases we remove the mean before creating the coeffs
%   illuminant structure
%     .wave  are wavelengths in nanometers
%     .data  are illuminant as a function of wavelength in energy units
%
% The full path to the output file is returned.
%
% The SPD of the data can be derived from the coefficients and basis
% functions using: 
%
%    spd = imageLinearTransform(mcCOEF,basis');
%
% See also: mcCreateMultispectralBases, CombineExposureColor
% See also: multicapture repository on github
%
% Copyright ImagEval Consultants, LLC, 2005.


if notDefined('mcCOEF'),  error('Coefficients required');     end
if notDefined('basis'),   error('Basis function required.');  end
if notDefined('comment'), comment = sprintf('Date: %s\n',date); end
if notDefined('illuminant'), error('Illuminant required'); end
if notDefined('fullName')
    fullName = vcSelectDataFile('stayput', ...
                        'w','mat','Save multispectral data file.');
end
 
% Write out the matlab data file with all of the key information needed.
% Sometimes we save out data approximated using only the SVD
% Other times, we use a principal component method and have an image mean
%
if notDefined('imgMean'), 
    save(fullName, 'mcCOEF', 'basis', 'comment', 'illuminant');
else
    save(fullName, 'mcCOEF', 'basis', 'imgMean', 'comment', 'illuminant');
end

end


