function RGB = imageSPD(SPD,wList,gam,row,col,displayFlag,xcoords,ycoords)
% Derive an RGB image from an SPD image
%
%     RGB = imageSPD(SPD,[wList],[gam],[row],[col],[displayFlag=1],[x coords],[y coords])
%
% The RGB image represents the appearance of the spectral power
% distribution (spd) data.  The RGB image is created by converting the
% image SPD to XYZ, and then converting the XYZ to sRGB format (xyz2srgb).
% The conversion to XYZ is managed so that the largest XYZ value is 1.
%
% The spd data should be in RGB format.  
%
% The routine can also be used to return the RGB values without displaying
% the data. 
%
% The image can be displayed with a spatial sampling grid overlaid to
% indicate the position on the optical image or sensor surface.
%
% wList: the sample wavelengths of the SPD 
%           (default depends on the number of wavelength samples)
% gam:   is the display gamma  (default = 1)
% row,col:  The image size (needed when the data are in XW format)
% displayFlag: <= 0 means no display. (default = 1, visible RGB and show)
%     Absolute value used to guide the type of rendering.  Thus,
%     if abs(displayFlag) = 
%              1:  Typical visible band to RGB
%              2:  
% x and y:  Coords (spatial positions) of the image points.
%   
% Examples:
%   imageSPD(spdData,[], 0.5,[],[],1);   % spdData is [r,c,w]
%   imageSPD(spdData,[], 0.6,128,128);   % spdData is [128*128,w]              
%   rgb = imageSPD(spdData,[], 0.5,[],[],0);  % rgb is calculated, but not displayed
%
% Copyright ImagEval Consultants, LLC, 2005.

if notDefined('gam'), gam = 1; end
if notDefined('displayFlag'), displayFlag = 1; end

if notDefined('wList')
    w = size(SPD, ndims(SPD));
    if     w == 31,  wList = (400:10:700); 
    elseif w == 301, wList = (400:1:700);
    elseif w == 37,  wList = (370:10:730);
    else 
        wList = sceneGet(vcGetObject('scene'),'wave');
        if length(wList) ~= w
            errordlg('Problem interpreting imageSPD wavelength list.');
            return;
        end
    end
end

% Convert the SPD data to a visible range image
if abs(displayFlag) == 1
    % RGB = imageSPD2RGB(SPD,wList,gam);
    XYZ = ieXYZFromPhotons(SPD,wList);
    XYZ = XYZ/max(XYZ(:));    
    RGB = xyz2srgb(XYZ);
elseif abs(displayFlag) == 2    % Gray scale image, used for SWIR, NIR
    RGB = zeros(row,col,3);
    RGB(:,:,1) = reshape(mean(SPD,3),row,col);
    RGB(:,:,2) = RGB(:,:,1);
    RGB(:,:,3) = RGB(:,:,1);
else
    error('Unknown display flag value: %d\n',displayFlag);
end

% Display the rendered RGB.  Sometimes we just return the RGB
% values. 
if displayFlag >= 1
    if ~isscalar(gam), RGB = RGB.^gam; end
    if notDefined('xcoords') || notDefined('ycoords')
        imagescRGB(RGB); axis image; 
    else
        RGB = RGB/max(RGB(:));
        RGB = ieClip(RGB,0,[]);
        imagesc(xcoords,ycoords,RGB);
        axis image;
        grid on; 
        set(gca,'xcolor',[.5 .5 .5]);
        set(gca,'ycolor',[.5 .5 .5]);
    end
end

end