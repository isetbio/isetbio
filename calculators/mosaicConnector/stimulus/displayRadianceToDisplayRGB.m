function [lRGBimage, sRGBimage, sRGBmaxContrast] = displayRadianceToDisplayRGB(radianceImage, spd)
    % Given a radiance image (M x N x W), and a set of display SPDs (W x 3),
    % compute the linear RGB primaries and sRGB values for that display
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    assert(size(spd,1) == wavelengthsNum, 'wavelengths do not match');
    
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
    linearRGB = spd \ radianceImage';
    minRGB = min(linearRGB(:));
    maxRGB = max(linearRGB(:));
    if (minRGB<0)
        linearRGB(linearRGB<0) = 0;
        minRGB = 0;
        fprintf('Stimulus not realizable. Will clip low end.\n');
    end
    if (maxRGB>1)
        linearRGB(linearRGB>1) = 1;
        maxRGB = 1;
        fprintf('Stimulus not realizable. Will clip high end.\n');
    end
    
    % Max contrast linear RGB
    linearRGBmaxContrast = (linearRGB-minRGB)/(maxRGB-minRGB);
    
    sRGB = lrgb2srgb(linearRGB');
    sRGBmaxContrast = lrgb2srgb(linearRGBmaxContrast');
    
    % Return in image format (MxNx3)
    lRGBimage = reshape(linearRGB, [rowsNum colsNum 3]);
    sRGBimage = reshape(sRGB, [rowsNum colsNum 3]);
    sRGBmaxContrast = reshape(sRGBmaxContrast, [rowsNum colsNum 3]);
end