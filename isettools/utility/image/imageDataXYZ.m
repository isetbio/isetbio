function dataXYZ = imageDataXYZ(vci, roiLocs)
% Return the XYZ values of display data 
%
% Syntax:
%   dataXYZ = imageDataXYZ(vci, [roiLocs])
%
% Description:
%    The linear primary data are contained in the result field. The XYZ
%    are computed by multiplying the monitor SPD (in energy) with these
%    linear RGB values of the display.
%
%    Typically, the entire image is returned in RGB (r, c, w)-format.
%
%    But, if roiLocs are passed in, the data are returned in XW format, 
%    with one X position for each roiLoc.
%
% Inputs:
%    vci     - 
%    roiLocs - (Optional). If present, return format is XW format
%
% Outputs:
%    dataXYZ - RGB or XW format of image (depending on roiLocs presence).
%
% Notes:
%    * [Note: JNM - This function doesn't work, as imageGet does not exist]
%    * [Note: JNM - This function is called in 0 places?]
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/06/17  jnm  Formatting

% Examples:
%{
    [val, vci] = vcGetSelectedObject('VCIMAGE');
    xyzRGB = imageDataXYZ(vci);
    roiLocs = vcROISelect(vci);
    xyzXW = imageDataXYZ(vci, roiLocs);
%}

if notDefined('roiLocs')
    % Get the rgb data. The result field contains linear RGB format for
    % the display.
    data = imageGet(vci, 'result');
        
    % Transform
    [data, r, c] = RGB2XWFormat(data);
    dataXYZ = imageRGB2XYZ(vci, data);
    dataXYZ = XW2RGBFormat(dataXYZ, r, c);
else
    % The data are returned in XW format    
    data = vcGetROIData(vci, roiLocs, 'result');
    dataXYZ = imageRGB2XYZ(vci, data);
end

end