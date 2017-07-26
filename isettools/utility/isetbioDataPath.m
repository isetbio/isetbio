function dataPath=isetbioDataPath()
%ISETBIODATAPATH  Return the path to the directory containing data bundled with isetbio.
%   dataPath=isetbioDataPath
%
%   See also: isetbioRootPath.

dataPath=fullfile(isetbioRootPath,'isettools','data',[]);

return
