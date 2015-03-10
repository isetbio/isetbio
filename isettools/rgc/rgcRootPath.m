function rootPath=rgcRootPath()
% Return the path to the root RGC directory
%
% This function must reside in the directory at the base of the ISET
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(rgcRootPath,'data')

rootPath=which('rgcRootPath');

[rootPath, ~, ~]=fileparts(rootPath);

return
