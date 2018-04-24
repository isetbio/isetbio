function RGB = imageSPD(SPD, wList, gam, row, col, displayFlag, ...
    xCoords, yCoords)
% Derive an RGB image from an SPD image
%
% Syntax:
%   RGB = imageSPD(SPD, [wList], [gam], [row], [col], [displayFlag], ...
%       [xCoords], [yCoords])
%
% Description:
%    The RGB image represents the appearance of the spectral power
%    distribution (spd) data. The RGB image is created by converting the
%    image SPD to XYZ, and then converting the XYZ to sRGB format
%    (xyz2srgb). The conversion to XYZ is managed so that the largest XYZ
%    value is 1.
%
%    The spd data should be in RGB format. 
%
%    The routine can also be used to return the RGB values without
%    displaying the data.
%
%    The image can be displayed with a spatial sampling grid overlaid to
%    indicate the position on the optical image or sensor surface.
%
%    Examples are located within the code. To access the examples, type
%    'edit imageSPD.m' into the Command Window.
%
% Inputs:
%    SPD         - The image in SPD format
%    wList       - (Optional) The sample wavelengths of the SPD. Default
%                  depends on the number of wavelength samples.
%    gam         - (Optional) The display gamma. Default is 1.
%    row         - (Optional) The image row size (needed for when data in
%                  XW). The default is pulled from provided image.
%    col         - (Optional) The image column size (needed for when data
%                  in XW). The default is pulled from the provided image.
%    displayFlag - (Optional) Whether or not the display is visible, and
%                  the type of rendering (color scheme). Default is 1,
%                  which means visible and RGB. The options are:
%                      0 - No Display (also true for values less than 0)
%                      1 - (Default) Typical visible band to RGB display
%                      2 - Gray-scale image. Used for SWIR and NIR.
%    xCoords     - (Optional) The x-axis spatial positions of image points.
%    yCoords     - (Optional) The y-axis spatial positions of image points.
%
% Outputs:
%    RGB         - The Resulting RGB image
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: JNM - The value 0 for displayFlag is clearly called out
%      above, yet the value was not present in the if statements for
%      calculations. I changed the if displayFlag == 1 to <= 1, so that for
%      a displayFlag of 0 it will be calculated as if RGB, but will not
%      display further on. This matches the language described above.]
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    12/08/17  jnm  Formatting & fix example
%    01/26/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    im = macbethChartCreate;
    spdData = im.data;

    % spdData is [r, c, w]
    imageSPD(spdData, [], 0.5, [], [], 1);

    % spdData is [128 * 128, w]
    imageSPD(spdData, [], 0.6, 128, 128);

    % rgb is calculated, but not displayed
    rgb = imageSPD(spdData, [], 0.5, [], [], 0);
%}

if notDefined('gam'), gam = 1; end
if notDefined('displayFlag'), displayFlag = 1; end
if notDefined('row'), [row, ~, ~] = size(SPD); end
if notDefined('col'), [~, col, ~] = size(SPD); end
if notDefined('wList')
    w = size(SPD, ndims(SPD));
    if  w == 31
        wList = (400:10:700); 
    elseif w == 301
        wList = (400:1:700);
    elseif w == 37
        wList = (370:10:730);
    else 
        wList = sceneGet(vcGetObject('scene'), 'wave');
        if length(wList) ~= w
            error('Problem interpreting imageSPD wavelength list.');
        end
    end
end

% Convert the SPD data to a visible color image (1) or gray scale (2) based
% on the absolute value.  The value is 0 or negative, so no imagesc()
% is called.
if abs(displayFlag) <= 1
    % Make RGB image, but do not display
    XYZ = ieXYZFromPhotons(SPD, wList);
    
    % We are considering getting rid of this normalization. The user
    % may want to set the relative intensity, so that two scenes with
    % different levels show up as lighter or darker RGB images as
    % well. By including this, we force all the images to be
    % normalized so that the brightest point is the same.
    XYZ = XYZ/max(XYZ(:));    
    RGB = xyz2srgb(XYZ);
    
elseif abs(displayFlag) == 2    
    % Gray scale image, used for SWIR, NIR
    RGB = zeros(row, col, 3);
    RGB(:, :, 1) = reshape(mean(SPD, 3), row, col);
    RGB(:, :, 2) = RGB(:, :, 1);
    RGB(:, :, 3) = RGB(:, :, 1);
else
    error('Unknown display flag value: %d\n', displayFlag);
end

% If the displayFlag is positive, then we want the data to be displayed,
% not just converted.  Hence the imagesc() calls
if displayFlag >= 1
    if ~isequal(gam, 1), RGB = RGB .^ gam; end
    
    % Choose the display method
    if notDefined('xcoords') || notDefined('ycoords')
        % No coordinates, just the image. This is used in the sceneWindow.
        imagescRGB(RGB);
        axis image; 
    else
        % This is used in pulldown callbacks and some scripts
        RGB = RGB / max(RGB(:));
        RGB = ieClip(RGB, 0, []);
        imagesc(xCoords, yCoords, RGB);
        axis image;
        grid on; 
        set(gca, 'xcolor', [.5 .5 .5]);
        set(gca, 'ycolor', [.5 .5 .5]);
    end
end

end