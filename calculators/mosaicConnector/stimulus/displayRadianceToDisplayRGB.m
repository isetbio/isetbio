function [lRGBimage, sRGBimage] = displayRadianceToDisplayRGB(radianceImage, spd)
    % Given a radiance image (M x N x W), and a set of display SPDs (W x 3),
    % compute the linear RGB primaries and sRGB values for that display
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    assert(size(spd,1) == wavelengthsNum, 'wavelengths do not match');
    
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
    linearRGB = spd \ radianceImage';
    if ((min(linearRGB(:))<0) || (max(linearRGB(:))>1))
        linearRGB(linearRGB<0) = 0;
        linearRGB(linearRGB>1) = 1;
        fprintf('Stimulus not realizable. Will clip.\n');
    end
    sRGB = lrgb2srgb(linearRGB');
    
    % Return in image format (MxNx3)
    lRGBimage = reshape(sRGB, [rowsNum colsNum 3]);
    sRGBimage = reshape(sRGB, [rowsNum colsNum 3]);
end