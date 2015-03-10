function fullName = sensorSaveImage(isa,fullName,dataType,gam,scaleMax)
%
%  fullpathName = sensorSaveImage(isa,fullName,dataType,gam,scaleMax);
%
% Save out an RGB image of the sensor image as a tiff file.  If the name
% is not passed in, then the user is queried to select the fullpath name
% of the output file.
% 
%Example:
%  sensorSaveImage(isa,fullpathname);
%
% See also: sensorData2Image, sensorGet(sensor,'rgb')
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('isa'), [~, isa] = vcGetSelectedObject('isa'); end
if notDefined('dataType'), dataType = 'volts'; end
if notDefined('gam'), gam = 1; end
if notDefined('scaleMax'), scaleMax = 1; end

if notDefined('fullName')
    fullName = vcSelectDataFile('session','w');
    if isempty(fullName), return;
    else
        [pathstr,name,~] = fileparts(fullName);
    end
    ext = '.tif';
    fullName = fullfile(pathstr,[name ext]);
end

% These are the displayed sensor data as an RGB image.
img = sensorData2Image(isa,dataType,gam,scaleMax);

% For internal display in Matlab, the img values range up to 255.  The
% imwrite requires doubles between [0,1].  Or uint8 between 0 and 255.
% We scale to 0,1.  This loses the absolute voltage level which must be
% recovered on the read side.
img = img/max(img(:));
imwrite(img,fullName,'tif');

end