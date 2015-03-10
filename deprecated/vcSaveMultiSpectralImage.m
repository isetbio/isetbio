function fullName = vcSaveMultiSpectralImage(imgDir,fname,mcCOEF,basis,basisLights,illuminant,comment,imgMean)
%
%   fullName = vcSaveMultiSpectralImage(imgDir,fname,mcCOEF,basis, ...
%                                   basisLights,illuminant,comment,imgMean)
%
%   Save a Matlab data file containing data for a multi-spectral image
%   
%   coefficients (RGB format), basis functions, illuminant information and
%   a comment. 
%
%   The full path to the data is returned in fullname.
%
%   The SPD of the data can be derived from the coefficients and basis
%   functions using: 
%
%   spd = rgbLinearTransform(mcCOEF,basis');
%

warning('Deprecated. Use ieSaveMultiSpectraImage instead');
if notDefined('basis'), error('Basis function required.');  end
if ~exist(imgDir, 'dir'), error('No such directory.'); end
if notDefined('comment')
    warning('Empty comment.');
    comment = sprintf('Date: %s\n',date);
end
if notDefined('basisLights')
    warning('No light description'); basisLights = []; 
end
if notDefined('illuminant')
    warning('No illuminant data'); illuminant = []; 
end

if isempty(fname)
    [fname, imgDir] = uiputfile('*-hdrs.mat', 'Enter HDRS file name');
    if isequal(fname,0) || isequal(imgDir,0)
        fullName = [];
        disp('HDRS file write canceled.')
        return;
    end
end
fullName = fullfile(imgDir,fname);

% Write out the matlab data file with all of the key information needed.
% Sometimes we save out data approximated usingly on the SVD
% Other times, we use a principal component method and have an image mean
if notDefined('imgMean')
    save fullName mcCOEF basis basisLights illuminant comment;
else
    save fullName mcCOEF basis basisLights illuminant imgMean comment;
end

end