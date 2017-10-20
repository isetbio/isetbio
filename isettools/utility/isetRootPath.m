function rootPath=isetRootPath()
%%isetRootPath  Return the path to the isetbio isettools directory
%
% Syntax:
%    rootPath = isetRootPath;
%
% Description:
%    Return base bath of the isetbio isettools directory.
%
% See also: isetbioRootPath, isetbioDataPath.

%% This is easy
rootPath = fullfile(isetbioRootPath,'isettools',[]);

return
