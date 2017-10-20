function dataPath=isetbioDataPath()
%%isetbioDataPath  Return the path to the oldsytle directory containing data bundled with isetbio.
%
% Syntax:
%    dataPath=isetbioDataPath;
%
% Description:
%    Return path to the isetbio data directory that is inside of isettools.
%
% Notes:
%    We are trying to make this directory go away, in favor for the data directory
%    and associated functions at the top level of isetbio.  But, this is a slow
%    process and we still need this for some older functions to work.
%
% See also: isetbioRootPath, isetRootPath.

dataPath=fullfile(isetbioRootPath,'isettools','data',[]);

return
