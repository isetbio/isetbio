function imT = imageFlip(im, flipType)
% Flip image data - updown or leftright 
%
% Syntax:
%   imT = imageFlip(im, flipType)
%
% Description:
%    The image data, im, size(im) = (r, c, w) can have any value for w. The
%    return, imT contains the data from im, but each color plane is flipped
%    according to the flipType: 'updown' or 'leftright'
%
% Inputs:
%    im       - The image data
%    flipType - (Optional) The flip type. Options are:
%               {'l', 'leftright'} - A horizontal flip. (Default)
%               {'u', 'updown'}    - A vertical flip. 
%
% Outputs:
%    imT      - The flipped image data.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    12/06/17  jnm  Formatting & fix example

% Examples:
%{
    list = 1:4;
    wave = 400:10:700;    % nanometers
    macbethReflectance = macbethReadReflectance(wave);
    XYZ = ieReadSpectra('XYZ',wave);
    lgt = ieReadSpectra('D65',wave);

    [U S V] = svd(macbethReflectance);
    W = S * V';
    mccXYZ = XYZ' * diag(lgt) * U(:,list) * W(list,:);

    imRGB = xyz2srgb(XW2RGBFormat(mccXYZ',4,6));
    imRGB = imageFlip(imRGB,'updown');

    imT1 = imageFlip(imRGB, 'upDown');
    imT2 = imageFlip(imRGB, 'leftRight');

    % subplot(2,3,1)
    % plot(imT1(:,:,1))
    % subplot(2,3,2)
    % plot(imT1(:,:,2))
    % subplot(2,3,3)
    % plot(imT1(:,:,3))
    % subplot(2,3,4)
    % plot(imT2(:,:,1))
    % subplot(2,3,5)
    % plot(imT2(:,:,2))
    % subplot(2,3,6)
    % plot(imT2(:,:,3))
%}

if ndims(im)~=3, error('Input must be rgb image (row x col x w)'); end
if notDefined('flipType'), flipType = 'l'; end

imT = zeros(size(im));
switch lower(flipType(1))
    case {'u', 'updown'}
        for ii=1:size(im, 3)
            imT(:, :, ii) = flipud(im(:, :, ii));
        end
        
    case {'l', 'leftright'}
        for ii=1:size(im, 3)
            imT(:, :, ii) = fliplr(im(:, :, ii));
        end
    otherwise
        error('unrecognized flip type');
end

end