 function [LMScontrastsImage, LMSexcitationsImage] = displayRadianceToLMS(radianceImage, coneFundamentals)
    rowsNum = size(radianceImage,1);
    colsNum = size(radianceImage,2);
    wavelengthsNum = size(radianceImage, 3);
    radianceImage = reshape(radianceImage, [rowsNum*colsNum wavelengthsNum]);
        
    LMSexcitations = radianceImage * coneFundamentals;
    backgroundLMS = LMSexcitations(1,:);
    
    LMScontrasts = bsxfun(@times, ...
        bsxfun(@minus, LMSexcitations, backgroundLMS), ...
        1./backgroundLMS);
    
    LMScontrastsImage = reshape(LMScontrasts, [rowsNum colsNum 3]);
    LMSexcitationsImage = reshape(LMSexcitations, [rowsNum colsNum 3]);
end
