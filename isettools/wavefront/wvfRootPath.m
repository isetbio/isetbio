function rootPath=wvfRootPath()
% Return the path to the root wavefront optics toolbox directory
%
% Syntax:
%   d = wvfRootPath
%
% Description:
%    This function must reside in the directory at the base of the
%    Wavefront Optics Tolbox directory structure. It is used to determine
%    the location of various sub-directories.
% 
% Inputs:
%    None.
%
% Outputs:
%    rootPath - The path to the root wavefront optics toolbox
%
% Notes:
%    * [Note: JNM - Is there a reason that the other arguments of fileparts
%      (fName, ext) are present, but not returned?]

% Examples:
%{
   d = wvfRootPath
%}

rootPath=which('wvfRootPath');

[rootPath, fName, ext]=fileparts(rootPath);

return