function dataPath = isetbioDataPath()
% Return path to oldstyle directory containing data bundled with isetbio
%
% Syntax:
%   dataPath = isetbioDataPath;
%
% Description:
%    Return path to the isetbio data directory that is inside of isettools.
%
% Inputs:
%    None.
%
% Outputs:
%    dataPath - The path to the old-style directory containing the data
%               bundled with isetbio.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - We are trying to make this directory go away, in favor
%      for the data directory and associated functions at the top level of
%      isetbio. But, this is a slow process and we still need this for some
%      older functions to work.]
%
% See Also:
%    isetbioRootPath, isetRootPath.

dataPath = fullfile(isetbioRootPath, 'isettools', 'data', []);

return
