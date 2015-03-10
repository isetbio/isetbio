function imageShowImage(vci,gam,trueSizeFlag,figNum)
% Display the image in the processor window.
%
%  imageShowImage(vci, [gam],[trueSizeFlag],[figNum])
%
% The RGB processed data in the virtual camera image (vci) are assumed to
% be linear RGB values.
%
% The processed RGB data are converted to sRGB.  The logic is that we have
% a display model for these data.  That display is stored in the vci. We
% convert the processed RGB data into XYZ values for that display.  We then
% convert the XYZ values into sRGB values, which we show to the user in the
% image processing window.
%
% Examples:
%   imageShowImage(vci{3},1/2.2)
%   imageShowImage(vci{3})
%
% See also:
%
% Copyright ImagEval Consultants, LLC, 2003.

% TODO
% I am concerned about the ordering of the ^gam and the scale
% operations.  Perhaps scaling should be first, and then the gamma.  As
% things stand, we apply gamma and then scale. 
% I am also unsured why the sRGB doesn't look better.  It looks good in the
% other windows. (BW).

if notDefined('vci'), cla; return;  end
if notDefined('gam'),  gam = 1; end
if notDefined('trueSizeFlag'), trueSizeFlag = 0; end
if notDefined('figNum'),  figNum = ieSessionGet('vcimageFigure'); end

% Bring up the figure
figure(figNum);

% Test and then convert the linear RGB values stored in result to XYZ
img = imageGet(vci,'result');
if isempty(img)
    cla; sprintf('There is no result image in vci.');
    return;
elseif max(img(:)) > imageGet(vci,'dataMax')
    error('Image max %.2f exceeds data max: %.2f\n',max(img(:)),dMax);
end

% Here is the key computational step.   We get the XYZ data from the vci.
% This calculation assumes that the processed RGB values are linear and the
% XYZ are calculated from the display primaries times these RGB.
% We then convert the XYZ to sRGB values.  We show those in the window.
img = imageGet(vci,'srgb');

if ismatrix(img),       vciType = 'monochrome';
elseif ndims(img) == 3, vciType = 'rgb';
else                    vciType = 'multisensor';
end

switch vciType
    case 'monochrome'
        colormap(gray(256));
        if gam ~= 1, imagesc(img.^(gam));
        else imagesc(img);
        end
    case 'rgb'
        if imageGet(vci,'scaleDisplay')
            % Use imagescRGB to render the RGB image.
            %  Prior to display negative values imagescRGB clips
            %  negative value and scales the result to a maximum of 1.
            if gam ~= 1, imagescRGB(img.^(gam));
            else         imagescRGB(img);
            end
            
        else
            % No scaling. There may be some negative numbers or numbers
            % > 1 because of processing noise and saturation + noise.
            img = ieClip(img,0,1);
            if gam ~= 1, image(img.^gam);
            else image(img);
            end
        end
    case 'multisensor'
        error('No display method for multisensor.');
end

axis image; axis off
if trueSizeFlag, truesize; end

end
