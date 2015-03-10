function rootPath=wvfRootPath()
% Return the path to the root wavefront optics toolbox directory
%
% This function must reside in the directory at the base of the Wavefront
% Optics Tolbox directory structure.  It is used to determine the location
% of various sub-directories.
% 
% Example:
%   d = wvfRootPath
%
% 

rootPath=which('wvfRootPath');

[rootPath,fName,ext]=fileparts(rootPath);

return