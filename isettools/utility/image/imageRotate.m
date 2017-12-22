function imT = imageRotate(im, rotType)
% Rotate image data - CW or CCW
%
% Syntax:
%   imT = imageRotate(im, rotType)
%
% Description:
%    The image data, im, size(im) = (r, c, w) can have any value for w.
%    imT contains the data from im, but each color plane is rotated.
%
%    We treat rotType = CW or CCW (ClockWise and CounterClockwise) as
%    special cases, using rot90.
%
%    If rotType is a number we call imrotate. Not sure this works properly
%    on all data. This feature is not yet used in ISET (I think).
%
% Inputs:
%    im      - 
%    rotType - (Optional). The rotation type. Default is 'ccw'. Options:
%           {'cw, 'clockwise'} - Clockwise rotation.
%           {'ccw, 'counterclockwise'} - Counter-clockwise rotation.
%
% Outputs:
%    imT     - The transposed image data
%
% Notes:
%    * [Note: JNM - Examples do not work with the image not instantiated. I
%      put in a macbethChart for the image, if this is not correct, please
%      let me know.]

% History:
%    xx/xx/09       Copyright ImagEval Consultants, LLC, 2009
%    12/07/17  jnm  Formatting

% Examples:
%{
    im = macbethChartCreate;
    imT1 = imageRotate(im.data, 'cw');
    imT2 = imageRotate(im.data, 'ccw');
    imT3 = imageRotate(im.data, 30);

    subplot(2,2,1)
    imagesc(im.data(:,:,25))
    subplot(2,2,2)
    imagesc(imT1(:,:,25))
    subplot(2,2,3)
    imagesc(imT2(:,:,25))
    subplot(2,2,4)
    imagesc(imT3(:,:,25))
%}

if ndims(im)~=3, error('Input must be rgb image (row x col x w)'); end
if notDefined('rotType'), rotType = 'ccw'; end

if isnumeric(rotType)
    tmp = imrotate(im(:, :, 1), rotType);
    [r, c, w] = size(tmp);
    clear tmp
    imT = zeros(r, c, w);
    for ii=1:size(im, 3)
        imT(:, :, ii) = imrotate(im(:, :, ii), rotType);
    end
else
    [r, c, w] = size(im);
    imT = zeros(c, r, w);
    switch lower(rotType)
        case {'cw', 'clockwise'}
            for ii=1:size(im, 3)
                imT(:, :, ii) = rot90(im(:, :, ii), -1);
            end
        case {'ccw', 'counterclockwise'}
            for ii=1:size(im, 3)
                imT(:, :, ii) = rot90(im(:, :, ii), 1);
            end
    end
end

end