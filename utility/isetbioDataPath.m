function dataPath = isetbioDataPath()
% Return path to data files bundled with ISETBio
%
% Syntax:
%   dataPath = isetbioDataPath;
%
% Description:
%    Return the canonical path to data files bundled with ISETBio.
%
% Inputs:
%    None.
%
% Outputs:
%    dataPath - The path to data files bundled with ISETBio.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    isetbioRootPath, isetRootPath.

dataPath = fullfile(isetbioRootPath, 'data', 'datafiles', []);

end
