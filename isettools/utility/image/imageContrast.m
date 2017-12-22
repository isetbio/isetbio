function cData = imageContrast(data)
% Compute the image contrast in an RGB style image
%
% Syntax:
%   cData = imageContrast(data)
%
% Description:
%    The image contrast is the intensity minus the mean divided by the
%    mean. Such an image always has zero mean, and can be useful for
%    computing the image MTF, discarding the DC term.
%
% Inputs:
%    data  - RGB style image data
%
% Outputs:
%    cData - Calculated image contrast data values
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/06/17  jnm  Formatting

cData = zeros(size(data));

for ii=1:size(data,3)
    m = mean(mean(data(:,:,ii)));
    cData(:,:,ii) = (data(:,:,ii) - m) / m;
end

end