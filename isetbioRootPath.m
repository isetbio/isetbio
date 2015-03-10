function rootPath=isetbioRootPath()
% Return the path to the root visual modeling directory
%
% This function must reside in the directory at the base of the vismodel
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(isetbioRootPath,'data')

rootPath=which('isetbioRootPath');

[rootPath,~,~]=fileparts(rootPath);

return
