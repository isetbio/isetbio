function imageType = ieImageType(fullName)
%Determine the type of image in a file
%
%   imageType = ieImageType(fullName)
%
% Check the file extensions to see if this is an RGB type
% file (e.g. jpeg,jpg,tif,gif,bmp).
%
% If not, read the directory name. If this contains one of the image type
% strings (see below), then return that string. Otherwise, we ask the user
% to identify the type of data. 
%
% Example:
%   fname =[isetRootPath,'\data\images\MultiSpectral\Fruit-hdrs.mat'];
%   ieImageType(fname)
%
% Copyright ImagEval Consultants, LLC, 2003.

[p,imageType,ext] = fileparts(fullName);
switch(lower(ext))
    case {'.jpg','.jpeg','.tif','.tiff','.bmp','.gif'}
        if strfind(fullName,'data\images\Targets'), 
            % Could be an EIA target.
            imageType = 'monochrome';
        elseif strfind(fullName,'data\images\Monochrome')
            imageType = 'monochrome';
        else       
            imageType = 'rgb';
        end
        return;
    otherwise
end

while ~isempty(imageType)
    imageType = lower(imageType);
    
    % If this is not a recognized type, ask the user.
    if (strcmp(imageType,'monochrome') || ...
            strcmp(imageType,'multispectral') || ...
            strcmp(imageType,'rgb') )
        return;
    end
    [p, imageType] = fileparts(p);
end

% If nothing in the path matches one of the known types, ask the user.
imageType = ieReadString('Enter file type {monochrome, rgb, or multispectral}','rgb');

end
